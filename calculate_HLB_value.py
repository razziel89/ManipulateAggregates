#this script is given an xyz-file and calculates the HLB-value using a VERY simple method

#import all the needed functions
from numpy import array, linalg, dot, arange, array_equal, inf, cross, sqrt
from numpy.linalg import norm
from math import pi, cos, sin, copysign, asin, acos
import sys
import re
import time

#this function returns the sign of a number by using lambda functions
sign = lambda x: copysign(1, x)

#this dictionary stores all masses of the atoms to be used in relative units
masses_dict={"H":1.0008, "O":16.0, "S":32.06, "N":14.01, "F":19.0, "Cl":35.45, "C":12.01, "Ni":58.69}
#this dictionary stores all weights to judge which atom is rather hydrophilic (high positive value) and which is not
weights_dict={"H":0, "O":1, "S":1, "N":1, "F":1, "Cl":1, "C":-0.4, "Ni":1}

from write_collection import *

#structure is a list of vectors
#normal_vector and coordinate define a plane (Hesse normal form)
#this function returns a list of indices of all atoms of the molecule that are one one side of the plane
#defined by normal_vector and coordinate
def part_molecule(structure,normal_vector,coordinate):
	return [i for i in range(0,len(structure)) if dot(array(structure[i])-array(coordinate),array(normal_vector)) >= 0 ]

#this function returns the index of the minimum element of a list liste
#name liste is German for list and has been chosen since list is a reserved word
def list_min(liste):
	return liste.index(min(liste))

#this function 
def pseudo_energy(structure,names,weights_dict,normal_vector,point_in_plane):
	hydrophilic_part_indices=part_molecule(structure,normal_vector,point_in_plane)
	lipophilic_part_indices=[i for i in range(0,len(structure)) if not i in hydrophilic_part_indices]
	return -sum([weights_dict[name] for name in [names[i] for i in hydrophilic_part_indices]])+sum([weights_dict[name] for name in [names[i] for i in lipophilic_part_indices]])

def define_new_z_axis(vector,new_axis,z_axis=array([0,0,1])):
        rotation_axis=cross(new_axis,z_axis)
        if norm(rotation_axis) == 0:
                return vector
        [x,y,z]=rotation_axis/norm(rotation_axis)
        c=dot(new_axis,z_axis)/(linalg.norm(new_axis)*linalg.norm(z_axis)) #cosine of the angle
        s=norm(cross(new_axis,z_axis))/(linalg.norm(new_axis)*linalg.norm(z_axis)) #sine of the angle
        t=1-c #definition that makes the calculation much faster
        R=array([[c+x**2*t,x*y*t-z*s,x*z*t+y*s],[x*y*t+z*s,c+y**2*t,y*z*t-x*s],[x*z*t-y*s,y*z*t+x*s,c+z**2*t]])
        return dot(R,vector)

#this function calculates the actual HLB-value
#this is done by moving a plane through the molecule that can go through any bond in any angle
#the standard values for the angles in which the plane can go through the bonds can be overwritten by
#declaring different phis and thetas in the function call (spherical co-ordinates are being used)
def HLB(structure,names,masses_dict,weights_dict,phis=arange(-pi,pi,2*pi/72),thetas=arange(0,pi,2*pi/72),relative_angles=False):
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
		pseudo_energy_list=[[pseudo_energy(structure,names,weights_dict,normal_vector,point_in_plane),normal_vector,point_in_plane] for point_in_plane in points_in_plane for normal_vector in normal_vectors]
	#find the index of the entry that contains the minimum energy
	min_index=list_min([liste[0] for liste in pseudo_energy_list])
	#extract the plane (as normal_vector and point_in_plane) that corresponds to the minimum in energy
	normal_vector=[liste[1] for liste in pseudo_energy_list][min_index]
	point_in_plane=[liste[2] for liste in pseudo_energy_list][min_index]
	#get the indices of the atoms in the hydrophilic part
	hydrophilic_part_indices=part_molecule(structure,normal_vector,point_in_plane)
	#calculate the HLB-value after Grffin's formula as 20 times the ratio of the mass in the hydrophilic part to the total mass of the molecule
	hlb_value=20*sum([masses_dict[name] for name in [names[i] for i in hydrophilic_part_indices]])/(sum([masses_dict[name] for name in names]))
	return hlb_value,normal_vector,point_in_plane
	

#try to open the given file for reading and throw error if not successfull
f=open(sys.argv[1],'r')

#read the lines in the given file into the variable lines
#and remove the trailing newline characters by using .rstrip()
lines = array([line.rstrip() for line in f])
#close the file descriptor
f.close()

#try to get the number of atoms in the molecule
#if this does not succeed, the file is probably not a valid
#xyz-file
try:
	nr_atoms=int(lines[0])
except ValueError:
	raise ValueError("this is probably not a valid xyz-file. The first line does not contain an integer.")

#the first two lines of an xyz file are not necessary, hence, they are removed
lines=array(lines[2:])

#for every line take the last three coloumns as co-ordinates for the atom 
#line.split() yields a list of the coloumns in line (whitespace separation)
#of those, take only the last ones line.split()[1:]
#map(float,L) gives a list whose elements are those of the list L but converted to floats
coordinates=array([map(float,line.split()[1:]) for line in lines])
#for every line take the element in the first coloumn as the name of the atom's element
names=[line.split()[0] for line in lines]

#calculate the HLB-value
#the plane that separates the hydrophilic and lipophilic parts is returned as well
anglestep=2*pi/72
(hlb_value,normal_vector,point_in_plane)=HLB(coordinates,names,masses_dict,weights_dict,thetas=array(list(arange(0,anglestep/2,anglestep))+list(arange(pi,pi-anglestep/2,anglestep))),relative_angles=True)
print hlb_value

#get the indices of the atoms in the hydrophilic part
hydrophilic_part_indices=part_molecule(coordinates,normal_vector,point_in_plane)
#get the lipophilic part as sonprising all atoms that are not in the hydrophilic part
lipophilic_part_indices=[i for i in range(0,len(coordinates)) if not i in hydrophilic_part_indices]
#print both parts
#'optimum.' is removed from the file name if present (done using regular expressions)
#[names[i] for i in hydrophilic_part_indices] gives a list of all element names of those atoms whose indices are in hydrophilic_part_indices
if len(re.findall(r"(?<=optimum\.).*",sys.argv[1])) != 0:
	filename=re.findall(r"(?<=optimum\.).*",sys.argv[1])[0]
else:
	filename=sys.argv[1]
print_xyz("hydrophilic_"+filename,[names[i] for i in hydrophilic_part_indices],[coordinates[i] for i in hydrophilic_part_indices])
print_xyz("lipophilic_"+filename,[names[i] for i in lipophilic_part_indices],[coordinates[i] for i in lipophilic_part_indices])
