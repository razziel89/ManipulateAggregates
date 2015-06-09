#-*- coding: utf-8 -*-
"""
Calculate the HLB-value using a simple method. Scan the molecule with a plane
that is allowed to cut midway through bonds with a maximum angle of 30° to
the bond.
"""
#import all the needed functions
import numpy as np
from math import pi, cos, sin, copysign, asin, acos

def pseudo_energy(point,vector,centers,potential):
    max_potential=np.max(potential)
    left_side=[True if np.dot(c-point,vector)>0 else False for c in centers]
    #left_energy = -sum((p for left,p in zip(left_side,potential) if left))
    #right_energy = -sum((max_potential-p for left,p in zip(left_side,potential) if not left))
    #print left_energy, right_energy
    #return left_energy+right_energy
    left_energy = -np.sum((p for left,p in zip(left_side,potential) if left))
    right_energy = np.sum((p-max_potential for left,p in zip(left_side,potential) if not left))
    #print left_energy, right_energy
    return left_energy+right_energy

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

def get_HLB(mol, nr_refinements):
    """
    """
    centers, potential = mol.get_vdw_surface_potential(nr_refinements=nr_refinements, weigh_by_surfacearea=True)
    centers=np.array(centers)
    abs_potential=np.fabs(np.array(potential))
    print abs_potential
    coordinates = np.array(mol.get_coordinates())
    bondmap = mol.get_bond_map(no_hydrogen=True)
    #first entry is the center of the bond
    #second entry is the direction of the bond
    bonds = (((coordinates[i]+coordinates[j])*0.5,coordinates[i]-coordinates[j]) for i,j in bondmap)

    #declare theta angles
    #TODO: rework this! This is not what it's supposed to be
    anglestep=2*pi/72
    thetas=np.concatenate((np.arange(-5*pi/180,5*pi/180,anglestep),np.arange(pi+5*pi/180,pi-5*pi/180,-anglestep)))
    phis=np.arange(-pi,pi,anglestep)

    #generate normal vectors allowed for the planes as vectors distributed over the border of a unit sphere
    normal_vectors=[(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) for theta in thetas for phi in phis]
    normal_vectors=[np.array(v) for v in set(normal_vectors)]

    #declare variables that will hold the optimum values
    min_energy = float("inf")
    normal_vector = None
    point_in_plane = None

    test_generator = ((point,define_new_z_axis(normal_vectors,direction)) for point,direction in bonds)

    for p in abs_potential:
        print p

    for point,vectors in test_generator:
        for vector in vectors:
            energy=pseudo_energy(point,vector,centers,abs_potential)
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
