"""
This script is able to arrange molecules arbitrarily.
Instead of doing so itself it rather provides the functions needed to do so.
Everything is wrapped in a class and uses openbabel data structures.
"""

import pybel as p
import sys
import numpy as np
import alphashapes as alpha

from electric_potential import potential_at_points

op=p.ob

def double_array(mylist):
    """Create a C array of doubles from a list."""
    c = op.doubleArray(len(mylist))
    for i,v in enumerate(mylist):
        c[i] = v
    return c

def drange(start, stop, step): #float range in steps of step
    """
    Like the built-in range but for floats.
    BEWARE: stop is included if hit!
    """
     r = start
     while r <= stop:
        yield r*step
        r += 1

def read_from_file(filename,fileformat="xyz"):
    """
    Read data from an xyz-file to a openbabel data structure.
    """
	return molecule(p.readfile(fileformat,filename).next().OBMol)

class molecule():
    """
    The molecule class aggregating all methods and data needed to describe and manipulate
    a molecule.
    """
	def __init__(self,mol,vector=None,axis=None,angle=None,ff='mmff94'):
        """
        Constructor.
        mol: the openbabel molecule data
        vector: displace the geometry that's read in from mol by this much (3-element vector)
        axis: rotate the geometry that's read in from mol around this axis (3-element vector)
        angle: the angle for the rotation
        ff: declare the forcefield associated with the molecule. Needed to get the energy and
            perform a simple forcefield geometry optimization.
        """
		self.mol=op.OBMol(mol)
#		op.OBGenericData.Clone(self.mol.OBGenericData,mol)
#		self.mol.CloneData(mol.GetData(op.OBGenericData_swigregister))
		if not ff in p.forcefields:
			return None
		self.ff=op.OBForceField.FindForceField(ff)
		self.ffname=ff
		if	self.ff.Setup(self.mol) == 0:
			return None
		#self.mol.Center()
		if not (axis == None or angle == None):
			self.rotate(axis,angle)
		if not vector == None:
			self.translate(vector)

	def duplicate(self):
        """
        Return a deep-copy of myself.
        """
		return molecule(self.mol,ff=self.ffname)

	def get_energy(self):
        """
        Get the energy associated with the current geometry for the current forcefield.
        """
		self.ff.Setup(self.mol)
		return self.ff.Energy()

	def optimize(self,steps=500):
        """
        Perform a sinple geometry optimization using the current forcefield.
        steps: number of optimization steps
        """
		p_tempmol=p.Molecule(self.mol)
		p_tempmol.localopt(forcefield=self.ffname,steps=steps)

	def set_bondlength(self,idx1,idx2,length,fix=None):
        """
        Adjust the length of a bond. If the bond connects two parts of a
        molecule that are otherwise not connected, those parts are moved
        with the respective atom. Otherwise, move only the 2 given atoms.
        There does not actually have to be a bond between the given atoms.

        Example:
        Original geometry:
         _   _
        |_|-|_|

        Adjust the middle bond to 3 times it's original length:
         _     _
        |_|---|_|

        Please note how both squares were moved along with the atoms
        comprising the bond.

        idx1: number of first atom that defines the bond
        idx2: number of second atom that defines the bond
        length: set the bond length to this value (probably in Angstrom)
        fix: if 1 or 2: keep the first or second atom fixed and move only the other.
             if None: move both by half the required distance
        """
		bond=op.OBBond()
		bond=self.mol.GetBond(idx1,idx2)
		if fix is None:
			bond.SetLength(length)
		elif fix == 1:
			bond.SetLength(self.mol.GetAtom(idx1),length)
		elif fix == 2:
			bond.SetLength(self.mol.GetAtom(idx2),length)

	def get_bondlength(self,idx1,idx2):
        """
        Get the length of a bond. There does not actually have to be a bond
        between the given atoms.

        idx1: number of first atom that defines the bond
        idx2: number of second atom that defines the bond
        """
		a1=op.OBAtom()
		a2=op.OBAtom()
		a1=self.mol.GetAtom(idx1)
		a2=self.mol.GetAtom(idx2)
		pos1=[a1.GetX(),a1.GetY(),a1.GetZ()]
		pos2=[a2.GetX(),a2.GetY(),a2.GetZ()]
		return ((a1.GetX()-a2.GetX())**2+(a1.GetY()-a2.GetY())**2+(a1.GetZ()-a2.GetZ())**2)**0.5
	
	def set_angle(self,idx1,idx2,idx3,angle):
        """
        Set the bond angle in rad. If the angle connects two parts of a
        molecule that are otherwise not connected, those parts are moved
        with the respective atom. See "set_bond" for a graphical example.

        idx1: number of first atom that defines the angle
        idx2: number of second atom that defines the angle
        idx3: number of third atom that defines the angle
        angle: set the angle to this value in rad
        """
                self.mol.SetAngle(idx1,idx2,idx3,angle)

	def get_angle(self,idx1,idx2,idx3):
        """
        Get the bond angle in rad. 

        idx1: number of first atom that defines the angle
        idx2: number of second atom that defines the angle
        idx3: number of third atom that defines the angle
        """
                return self.mol.GetAngle(idx1,idx2,idx3)

	def set_dihedral(self,idx1,idx2,idx3,idx4,angle):
        """
        Set the dihedral angle in rad. If the angle connects two parts of a
        molecule that are otherwise not connected, those parts are moved
        with the respective atom. See "set_bond" for a graphical example.

        idx1: number of first atom that defines the dihedral angle
        idx2: number of second atom that defines the dihedral angle
        idx3: number of third atom that defines the dihedral angle
        idx4: number of fourth atom that defines the dihedral angle
        angle: set the angle to this value in rad
        """
		self.mol.SetDihedralAngle(idx1,idx2,idx3,idx4,angle)

	def get_dihedral(self,idx1,idx2,idx3,idx4):
        """
        Get the dihedral angle in rad. 

        idx1: number of first atom that defines the dihedral angle
        idx2: number of second atom that defines the dihedral angle
        idx3: number of third atom that defines the dihedral angle
        idx4: number of fourth atom that defines the dihedral angle
        """
		return self.mol.GetDihedralAngle(idx1,idx2,idx3,idx4)

	def rotate(self,axis,angle):
        """
        Rotate the molecule around an axis by an angle.

        axis: rotate the geometry around this axis (3-element vector)
        angle: the angle for the rotation
        """
		matrix=op.matrix3x3();
		matrix.RotAboutAxisByAngle(op.vector3(double_array(axis)),angle)
		array=double_array([0]*9)
		matrix.GetArray(array)
		self.mol.Rotate(array)

	def rotate_main(self,axis_index,angle):
        """
        Rotate the molecule around one of its main axis by an angle.

        axis_index: 1, 2 or 3: the index of the main axis to rotate around
        angle: the angle for the rotation
        """
		self.mol.Rotate(axis_index,angle)

	def vdw_check(self):
        """
        Check whether any two atoms are closer together than the sum of
        their van-der-Vaals radii. Perform this check only for atoms that
        are not connected by an arbitrary number of bonds. Hence, this
        only makes sense for aggregates.
        """
		return self.mol.IsGoodVDW()

	def translate(self,vector):
        """
        Translate the molecule in a given direction.

        vector: 3-element vector that is added to every atom's coordinate
        """
		self.mol.Translate(double_array(vector))

	def append(self,mol,vector=[0,0,0],axis=[1,0,0],angle=0):
        """
        Append a molecule to the current one. Before appending, translate and rotate
        the part that is to be appended.

        mol: molecule to be appended (of type molecule)
        vector: 3-element vector that is added to the coordinate of
                every atom that will be appended
        axis: rotate the to be appended geometry around this axis (3-element vector)
        angle: the angle for the rotation
        """
		self.mol.AppendMolecule(mol.mol,double_array(vector),double_array(axis),angle)

	def write(self,filename,overwrite='False',fileformat='xyz'):
        """
        Write the data of the molecule to disk. 

        filename: name of the file INCLUDING the extension
        overwrite: shall the output file be overwritten or not
        fileformat: output file format (anything that openbabel can write)
        """
		p.Molecule(self.mol).write(fileformat,filename,overwrite=overwrite)

	def align(self,point,main1,main2):
        """
        Align the first two main axes of a molecule to the two given axis and
        move the center to the given coordinate.

        point: 3-element sequence defining the new center of the molecule 
               (not mass weighed)
        main1: 3-element sequence defining the new 1st main axis
        main2: 3-element sequence defining the new 2nd main axis

        """
		if (sum(main1)>0 and sum(main2)>0):
			self.mol.Align(double_array(point),double_array(main1),double_array(main2))

	def part_molecule(self,normal_vector,coordinate):
        """
        Get a molecule containing all those atoms that are one one side of a plane given
        in the Hessian normal form.

        normal_vector: normal vector of Hessian normal form (3-element list)
        coordinate: 3d-Cartesian coordinates of one point in the plane
        """
		tempmol=op.OBMol()
		self.mol.PartMolecule(tempmol,double_array(normal_vector),double_array(coordinate))
		return tempmol

	def write_part(self,filename,normal_vector,coordinate,side='left',overwrite='False',fileformat='xyz'):
        """
        Write the data of the molecule to disk. Write only those atoms to disk that are
        one one side of a plane given in the Hessian normal form.

        filename: name of the file INCLUDING the extension
        normal_vector: normal vector of Hessian normal form (3-element list)
        coordinate: 3d-Cartesian coordinates of one point in the plane
        side: If 'left', print those atoms on the side where the normal vector points.
              If anything else, print all those on the opposite side.
        overwrite: shall the output file be overwritten or not
        fileformat: output file format (anything that openbabel can write)
        """
        if side=='left':
            p.Molecule(self.part_molecule(normal_vector,coordinate)).write(fileformat,filename,overwrite=overwrite)
        else:
            normal_vector=[-i for i in normal_vector]
            p.Molecule(self.part_molecule(normal_vector,coordinate)).write(fileformat,filename,overwrite=overwrite)
	
#	def HLB_value(self, nr_points=100):
#            a=op.OBAtom()
#            partialcharges=[0.0]*self.mol.NumAtoms()
#            coordinates=[None]*self.mol.NumAtoms()
#            vdw_radii=[0.0]*self.mol.NumAtoms()
#            masses=[0.0]*self.mol.NumAtoms()
#            non_hydrogen_index=[]
#            for idx in range(1,self.mol.NumAtoms()+1):
#                a = self.mol.GetAtom(idx)
#                partialcharges[idx-1] = a.GetPartialCharge()
#                coordinates[idx-1] = [a.GetX(),a.GetY(),a.GetZ()]
#                vdw_radii[idx-1] = op.etab.GetVdwRad(a.GetAtomicNum())
#                masses[idx-1] = op.etab.GetMass(a.GetAtomicNum())
#            points = [ point for i in range(0,self.mol.NumAtoms()) for point in alpha.sphere_distribution(nr_points, i, coordinates, vdw_radii) ]
#            [trig_centres, areas] = alpha.triangulated_alphashape(points)
#            potential = potential_at_points(trig_centres, partialcharges, coordinates)
