#-*- coding: utf-8 -*-
"""
Calculate the HLB-value using a simple method. Scan the molecule with a plane
that is allowed to cut midway through bonds with a maximum angle to
the bond.
"""
#import all the needed functions
import numpy as np
from math import pi, cos, sin, copysign, asin, acos

class NoBondsException(Exception):
    pass

def pseudo_energy(max_potential,point,vector,centers,potential,weights):
    """
    Estimate the energy the system would have if separated by a plane on whose
    left (IN the direction of the normal vector) there is an unpolar solvend
    and on whose right there is a polar solvent, like H2O. That means,
    aggregating potential!=0 on the right side is good whereas aggregating
    potential close to 0 is good on the left.

    max_potential: the maximum value of the potential "potential"
    point: an arbitrary point on the plane that separates the system
    vector: the normal vector of the plane that separates the system.
            This need not be normalized but MUST HAVE A NON-VANISHING NORM.
    centers: Cartesian coordinates where the potentials are located.
             Must be a 1D-numpy array of vectors.
    potential: a numpy array containing the potential values.
               Has to be 1D and contain only positive values.
    weights: if some potentials shall count more than others, for instance
             one bit should in every case be on a specific side, than declare
             something other than None for this.
             Here, since each potential is in the center of a triangle, this is
             the area of the triangle.
    """
    left_side=np.dot(centers-point,vector)>0
    right_side=np.invert(left_side)
    energy_left = np.average(potential,weights=weights*left_side)
    energy_right=-np.average(potential,weights=weights*right_side)
    return energy_left+energy_right

def define_new_z_axis(vector_list,new_axis,z_axis=np.array([0,0,1])):
    """
    Compute the rotation matrix that can rotate the current z-axis to a new axis
    and apply this rotation matrix to a list of vectors.

    vector_list: vectors that shall be rotated
    new_axis: new z-axis
    z_axis: current z-axis
    """
    rotation_axis=np.cross(new_axis,z_axis)
    if np.linalg.norm(rotation_axis) == 0:
        return vector
    [x,y,z]=rotation_axis/np.linalg.norm(rotation_axis)
    c=np.dot(new_axis,z_axis)/(np.linalg.norm(new_axis)*np.linalg.norm(z_axis)) #cosine of the angle
    s=np.linalg.norm(np.cross(new_axis,z_axis))/(np.linalg.norm(new_axis)*np.linalg.norm(z_axis)) #sine of the angle
    t=1-c #definition that makes the calculation much faster
    R=np.array([[c+x**2*t,x*y*t-z*s,x*z*t+y*s],[x*y*t+z*s,c+y**2*t,y*z*t-x*s],[x*z*t-y*s,y*z*t+x*s,c+z**2*t]])
    return (np.dot(R,vector) for vector in vector_list)

def get_HLB(mol, nr_refinements,anglestep_degree=5,anglerange_degree=10,no_hydrogen=True):
    """
    Compute the HLB-value of a molecule of the molecule class defined in manipulate_molecules module.
    Returns the hlb value, the normal vector and a point in the plane that divides the molecule.

    mol: molecule of the aforementioned class
    nr_refinements: how many refinement steps to use during skin surface generation
    anglestep_degree: in what steps shall the plane be allowed to vary
    anglerange_degree: plus-minus this many degrees is the normal vector allowed to
                       stray from the bond vector
    no_hydrogen: whether or not to irgnore bonds to hydrogens as cutable bonds
    """
    surfaceareas=[]
    centers, potential = mol.get_vdw_surface_potential(nr_refinements=nr_refinements, weights=surfaceareas)
    surfaceareas=np.array(surfaceareas)
    centers=np.array(centers)
    abs_potential=np.fabs(np.array(potential))
    coordinates = np.array(mol.get_coordinates())
    bondmap = mol.get_bond_map(no_hydrogen=no_hydrogen)
    if len(bondmap)==0:
        raise NoBondsException("Number of bonds through which cutting would be allowed is 0!")

    #first entry is the center of the bond
    #second entry is the direction of the bond
    bonds = (((coordinates[i]+coordinates[j])*0.5,coordinates[i]-coordinates[j]) for i,j in bondmap)

    #declare theta angles
    #TODO: rework this! This is not what it's supposed to be
    anglestep=anglestep_degree*pi/180
    anglerange=anglerange_degree*pi/180
    thetas=np.concatenate((np.arange(-anglerange,anglerange,anglestep),np.arange(pi+anglerange,pi-anglerange,-anglestep)))
    phis=np.arange(-pi,pi,anglestep)

    #generate normal vectors allowed for the planes as vectors distributed over the border of a unit sphere
    normal_vectors=[(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) for theta in thetas for phi in phis]
    normal_vectors=[np.array(v) for v in set(normal_vectors)]

    #declare variables that will hold the optimum values
    min_energy = float("inf")
    normal_vector = None
    point_in_plane = None

    test_generator = ((point,define_new_z_axis(normal_vectors,direction)) for point,direction in bonds)

    max_potential=np.max(potential)
    
    for point,vectors in test_generator:
        for vector in vectors:
            energy=pseudo_energy(max_potential,point,vector,centers,abs_potential,surfaceareas)
            if energy<min_energy:
                min_energy=energy
                normal_vector=vector
                point_in_plane=point
    #            print energy,normal_vector,point_in_plane
    
    lipomass=sum(mol.part_molecule_mol(normal_vector,point_in_plane).get_masses())
    totalmass=sum(mol.get_masses())
    #calculate the HLB-value after Grffin's formula as 20 times the ratio of the mass in the hydrophilic part to the total mass of the molecule
    hlb_value=20*(1-lipomass/totalmass)

    #part molecule
    #[i for i in range(0,len(structure)) if np.dot(np.array(structure[i])-np.array(coordinate),np.array(normal_vector)) >= 0 ]

    return hlb_value,normal_vector,point_in_plane
