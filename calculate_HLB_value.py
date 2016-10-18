#-*- coding: utf-8 -*-
"""
Calculate the HLB-value using a simple method. Scan the molecule with a plane
that is allowed to cut midway through bonds with a maximum angle of 30° to
the bond. Assign positive values to all hetero atoms and negative values to C and H. Try to maximize the number of heteroatoms on one side and the number
of C and H on the other side.
"""
#This file is part of ManipulateAggregates.
#
#Copyright (C) 2016 by Torsten Sachse
#
#ManipulateAggregates is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ManipulateAggregates is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.

from numpy import array, linalg, dot, arange, array_equal, inf, cross, sqrt
from numpy.linalg import norm
from math import pi, cos, sin, copysign, asin, acos
import sys
import re
import time

from collection.read import *
from collection.write import *

global masses_dict, weights_dict
#this dictionary stores all masses of the atoms to be used in relative units
masses_dict={"H":1.0008, "O":16.0, "S":32.06, "N":14.01, "F":19.0, "Cl":35.45, "C":12.01, "Ni":58.69}
#this dictionary stores all weights to judge which atom is rather hydrophilic (high positive value) and which is not
weights_dict={"H":0, "O":1, "S":1, "N":1, "F":1, "Cl":1, "C":-0.4, "Ni":1}

#this function returns the sign of a number by using lambda functions
sign = lambda x: copysign(1, x)

def part_molecule(structure,normal_vector,coordinate):
    """
    this function returns a list of indices of all atoms of the molecule that are one one side of the plane
    defined by normal_vector and coordinate

    structure: list of vectors
    normal_vector: normal vector for Hessian normal form
    coordinate: point in the plane for Hessian normal form
    """
    return [i for i in range(0,len(structure)) if dot(array(structure[i])-array(coordinate),array(normal_vector)) >= 0 ]

def list_min(liste,key=lambda x:x):
    """
    this function returns the index of the minimum element of a list liste
    name liste is German for list and has been chosen since list is a reserved word
    """
    li=[key(l) for l in liste]
    return li.index(min(li))

def pseudo_energy(structure,names,weights_dict,normal_vector,point_in_plane):
    """
    With the given structure, calculate the "energy" associated with it if the
    plane, given in its Hessian normal form, cuts through the molecule. The value
    is high and positive, if atoms with a positive weight are on the side where
    the normal vector points and atoms with negative weights are on the opposite
    side of the plane.

    structure: a list of 3-element vectors containing the cartesian coordinates
               of all the atoms in the structure
    names: a list containing the element names of all the atoms in the structure
    weights_dict: a dictionary that has a weight assigned to each element name
    normal_vector: the normal vector of the plane in its Hessian normal form
    point_in_plane: the point (3 cartesian coordinates) that lies in the plane
    """
    hydrophilic_part_indices=part_molecule(structure,normal_vector,point_in_plane)
    lipophilic_part_indices=part_molecule(structure,-array(normal_vector),point_in_plane)
    return -sum([weights_dict[name] for name in [names[i] for i in hydrophilic_part_indices]])+sum([weights_dict[name] for name in [names[i] for i in lipophilic_part_indices]])

def define_new_z_axis(vector,new_axis,z_axis=array([0,0,1])):
    """
    Compute the rotation matrix that can rotate the current z-axis to a new axis
    and apply this rotation matrix to a vector.

    vector: the vector that shall be rotated
    new_axis: new z-axis
    z_axis: current z-axis

    Example:
    Old state (z: old z-axis, n: new z-axis, vb: vector before, va: vector after):

    z   vb
    ↑◜◝↗
    |φ╱
    |╱
    |------→ n

    New state:

    z   
    ↑
    |
    |
    |------→ n
     ╲φ◝
      ╲◞
       ↘ va
        
    """
    rotation_axis=cross(new_axis,z_axis)
    if norm(rotation_axis) == 0:
            return vector
    [x,y,z]=rotation_axis/norm(rotation_axis)
    c=dot(new_axis,z_axis)/(linalg.norm(new_axis)*linalg.norm(z_axis)) #cosine of the angle
    s=norm(cross(new_axis,z_axis))/(linalg.norm(new_axis)*linalg.norm(z_axis)) #sine of the angle
    t=1-c #definition that makes the calculation much faster
    R=array([[c+x**2*t,x*y*t-z*s,x*z*t+y*s],[x*y*t+z*s,c+y**2*t,y*z*t-x*s],[x*z*t-y*s,y*z*t+x*s,c+z**2*t]])
    return dot(R,vector)

def HLB(structure,names,masses_dict,weights_dict,phis=arange(-pi,pi,2*pi/72),thetas=arange(0,pi,2*pi/72),relative_angles=False,pseudo_energy=pseudo_energy):
    """
    this function calculates the actual HLB-value
    this is done by moving a plane through the molecule that can go through any bond in any angle
    the standard values for the angles in which the plane can go through the bonds can be overwritten by
    declaring different phis and thetas in the function call (spherical co-ordinates are being used)

    structure: a list of 3-element vectors containing the cartesian coordinates
               of all the atoms in the structure
    names: a list containing the element names of all the atoms in the structure
    masses_dict: a dictionary that has a mass assigned to each element name in atomic units
    weights_dict: a dictionary that has a weight assigned to each element name
    phis: a list containing all the angles in the x-y-plane the plane should screen
          if relative_angles == true, this will correspond to the new x-y-plane after
          the rotation
    thetas: same as phis but for the theta angle of spherical coordinates
    relative_angles: if true, realign the z-axis to the bond vector prior to scanning
    """
    #remove from the list of coordinates those that belong to hydrogen atoms since
    #the plane is not allowed to go though a bond with hydrogen
    structure_no_H=[structure[i] for i in range(0,len(names)) if names[i]!="H"]

    #generate normal vectors allowed for the planes as vectors distributed over the border of a unit sphere
    normal_vectors=[array([sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]) for theta in thetas for phi in phis]
    normal_vectors=[array(v) for v in set([tuple(vector) for vector in normal_vectors])]

    #the plane has to cut somewhere through a bond
    #the cutting point is taken as the middle of a bond
    #a bond is considered to exist if the two atom sare closer together than 2.55 Angstroms
    #set removes duplicates but this cannot be done for lists of coordinates only for tuples
    #since lists cannot be hashed
    if relative_angles:
        #find out the coordinates of those atoms that belong to bonds and remove duplicates
        bond_atoms=[(tuple(a),tuple(b)) for a in structure_no_H for b in structure_no_H if norm(a-b) < 2.55]
        bond_atoms=[[array(coord) for coord in two_coords] for two_coords in set(bond_atoms)]
        #the function pseudo_energy returns an energy value per plane the lowest of which corresponds to the
        #optimum position of an interface between a hydrophilic and lipophilic solvent
        #this can be used to guess the way the molecule arranges itself at such an interface
        pseudo_energy_list=[]
        for bond in bond_atoms:
            bond_direction=bond[0]-bond[1]
            point=sum(bond)/2
            for normal_vector in normal_vectors:
                new_normal_vector=define_new_z_axis(normal_vector,bond_direction)
                pseudo_energy_list.append([pseudo_energy(structure,names,weights_dict,new_normal_vector,point),new_normal_vector,point])
    else:
        points_in_plane=[tuple((a+b)/2) for a in structure_no_H for b in structure_no_H if norm(a-b) < 2.55]
        #recast the coordinates to numpy arrays to allow easy processing
        points_in_plane=[array(coord) for coord in set(points_in_plane)]
        pseudo_energy_list=[]
        for point_in_plane in points_in_plane:
            for normal_vector in normal_vectors:
                pseudo_energy_list.append([pseudo_energy(structure,names,weights_dict,normal_vector,point_in_plane),normal_vector,point_in_plane])

    #find the index of the entry that contains the minimum energy
    min_index=list_min(pseudo_energy_list,key=lambda x:x[0])
    #extract the plane (as normal_vector and point_in_plane) that corresponds to the minimum in energy
    normal_vector=pseudo_energy_list[min_index][1]
    point_in_plane=pseudo_energy_list[min_index][2]
    #get the indices of the atoms in the hydrophilic part
    hydrophilic_part_indices=part_molecule(structure,normal_vector,point_in_plane)
    #calculate the HLB-value after Grffin's formula as 20 times the ratio of the mass in the hydrophilic part to the total mass of the molecule
    hlb_value=20*sum([masses_dict[name] for name in [names[i] for i in hydrophilic_part_indices]])/(sum([masses_dict[name] for name in names]))
    return hlb_value,normal_vector,point_in_plane

def main():
    global masses_dict, weights_dict

    #read in the data
    names,coordinates=read_xyz(sys.argv[1])
    
    #declare theta angles
    anglestep=2*pi/72
    thetas=array(list(arange(0,anglestep/2,anglestep))+list(arange(pi,pi-anglestep/2,anglestep)))

    #calculate the HLB-value
    #the plane that separates the hydrophilic and lipophilic parts is returned as well
    #angles are given relative to the bond
    (hlb_value,normal_vector,point_in_plane)=HLB(coordinates,names,masses_dict,weights_dict,thetas=thetas,relative_angles=True)
    print hlb_value
    
    #get the indices of the atoms in the hydrophilic part
    hydrophilic_part_indices=part_molecule(coordinates,normal_vector,point_in_plane)
    lipophilic_part_indices=part_molecule(coordinates,-array(normal_vector),point_in_plane)

    #print both parts
    #'optimum.' is removed from the file name if present (done using regular expressions)
    #[names[i] for i in hydrophilic_part_indices] gives a list of all element names of those atoms whose indices are in hydrophilic_part_indices
    if len(re.findall(r"(?<=optimum\.).*",sys.argv[1])) != 0:
        filename=re.findall(r"(?<=optimum\.).*",sys.argv[1])[0]
    else:
        filename=sys.argv[1]

    print_xyz("hydrophilic_"+filename,[names[i] for i in hydrophilic_part_indices],[coordinates[i] for i in hydrophilic_part_indices])
    print_xyz("lipophilic_"+filename,[names[i] for i in lipophilic_part_indices],[coordinates[i] for i in lipophilic_part_indices])

    
if __name__ == "__main__":
    main()
