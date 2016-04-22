"""
This script is able to arrange molecules arbitrarily.
Instead of doing so itself it rather provides the functions needed to do so.
Everything is wrapped in a class and uses openbabel data structures.
"""

global supported
supported={}

class ManipulateMoleculesError(Exception):
    pass

class MissingModuleError(ManipulateMoleculesError):
    pass

class FiletypeException(ManipulateMoleculesError):
    pass

class OpenBabelError(ManipulateMoleculesError):
    pass

import os
import re
import sys
import copy

try:
    import pybel as p
except ImportError as e:
    raise MissingModuleError("Pybel could not be imported. Please install openbabel with Python bindings.",e)
try:
    import numpy as np
    supported["numpy"]=(True,)
except ImportError as e:
    supported["numpy"]=(False,e)

try:
    import FireDeamon as fd
    supported["FireDeamon"]=(True,)
except ImportError as e:
    supported["FireDeamon"]=(False,e)

import electric_potential as ep

op=p.ob

#copy all formats known to openbabel into a new dictionary
filetypedict={entry:entry for entry in p.informats}
#add some custom entries
filetypedict["mop"]="mopin"

def _double_array(mylist):
    """Create a C array of doubles from a list."""
    c = op.doubleArray(len(mylist))
    for i,v in enumerate(mylist):
        c[i] = v
    return c

def guess_format(filename):
    """
    Try to guess the file format. If the filename contains no dot,
    try to match the whole filename agains the database. Useful
    for e.g. CONTCAR files that seldomnly contain a dot.
    """
    if re.match("^[^.]+$",filename)==None:
        match=re.search('(?<=\.).*$',filename)
        extension=match.string[match.start():match.end()]
    else:
        extension=filename
    try:
        filetype=filetypedict[extension]
    except KeyError as e:
        raise FiletypeException("Filetype of file "+filename+" not known to openbabel.",e)
    return filetype

def read_from_file(filename,fileformat=None,conf_nr=1,ff='mmff94'):
    """
    Read data from a file to a openbabel data structure. Guess filetype
    if none present.

    fileformat: guess type if None, otherwise use specified filetype
    conf_nr: can be a single number specifying the conformer to load
             or can be an iterable returning the indices of those
             conformers to load. Special keyword 'all' will return all
             conformers in file.
    """
    if fileformat==None:
        fileformat=guess_format(filename)
    if re.match("~.*/",filename):
        homedir=filename.split("/")[0]
        if re.match("^~$",homedir):
            filename=re.sub("^~",os.environ["HOME"],filename)
        else:
            username=homedir.split("~")[1]
            homedir="/".join(os.environ["HOME"].split("/")[:-1])
            filename=re.sub("^~",homedir+"/",filename)

    fileinfo = {'name':filename, 'format':fileformat, 'conf_nr':conf_nr, 'ff': ff}
    if conf_nr==1:
        mol = molecule(p.readfile(fileformat,filename).next().OBMol, ff=ff, fileinfo=fileinfo)
    else:
        try:
            conf_nr_iter = iter(conf_nr)
            iterable = True
        except TypeError:
            iterable = False
        conformers=[m for m in p.readfile(fileformat,filename)]
        if conf_nr=='all':
            conf_nr=range(1,len(conformers)+1)
            iterable=True
        if iterable:
            if max(conf_nr)>len(conformers):
                raise ValueError("You requested conformer number %d but there are only %d present in the file."%(max(conf_nr),len(conformers)))
            return [molecule(conformers[i-1].OBMol) for i in conf_nr]
        else:
            if conf_nr>len(conformers):
                raise ValueError("You requested conformer number %d but there are only %d present in the file."%(conf_nr,len(conformers)))
            mol = molecule(conformers[conf_nr-1].OBMol, ff=ff, fileinfo=fileinfo)
    return mol

def _RotMatrixAboutAxisByAngle(axis,angle):
    """
    Taken from openbabel from file matrix3x3.cpp, method RotAboutAxisByAngle.
    Generate a rotation matrix about an arbitrary axis by an arbitrary angle.
    Angle has to be in radians.
    """
    mat = np.identity(3,dtype=float)
    theta = angle;
    s = np.sin(theta);
    c = np.cos(theta);
    t = 1.0 - c;

    vtmp = np.array(axis)
    if not len(vtmp.shape)==1 and vtmp.shape[0]==3:
        raise ValueError("Given axis must have shape (3,) but it has shape "+str(vtmp.shape))
    if np.linalg.norm(vtmp)>0.001:
        vtmp /= np.linalg.norm(vtmp);

        x,y,z = vtmp

        mat[0][0] = t*x*x + c;
        mat[0][1] = t*x*y + s*z;
        mat[0][2] = t*x*z - s*y;

        mat[1][0] = t*x*y - s*z;
        mat[1][1] = t*y*y + c;
        mat[1][2] = t*y*z + s*x;
        
        mat[2][0] = t*x*z + s*y;
        mat[2][1] = t*y*z - s*x;
        mat[2][2] = t*z*z + c;

    return mat

def _VectorAngle(vector1,vector2):
    """
    Return the angle between two vectors in radians.
    """
    v1 = np.array(vector1)
    v2 = np.array(vector2)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dp = np.dot(v1,v2)
    if dp < -1.0:
        dp = -1.0
    elif dp > 1.0:
        dp = 1.0
    return np.arccos(dp)

class molecule():
    """
    The molecule class aggregating all methods and data needed to describe and manipulate
    a molecule.
    """
    def __init__(self,mol,vector=None,axis=None,angle=None,ff='mmff94',fileinfo={},charge_method='gasteiger'):
        """
        Constructor.
        mol: the openbabel molecule data
        vector: displace the geometry that's read in from mol by this much (3-element vector)
        axis: rotate the geometry that's read in from mol around this axis (3-element vector)
        angle: the angle for the rotation
        ff: declare the forcefield associated with the molecule. Needed to get the energy and
            perform a simple forcefield geometry optimization.
        fileinfo: a dictionary with keys 'name' and 'format' detailing the filename and format,
                  respectively. Also, the keys 'conf_nr' and 'ff' for the used conformer number
                  and force field are required.
        """
        self.mol=op.OBAggregate(mol)
        self.fileinfo=copy.deepcopy(fileinfo)
        self.charge_method=charge_method
        self.ffname=ff
        if not ff==None:
            if not ff in p.forcefields:
                print >> sys.stderr, "Force field not known to openbabel."
                #return None
            self.ff=op.OBForceField.FindForceField(ff)
            if  self.ff.Setup(self.mol) == 0:
                print >> sys.stderr, "Force field could not be set-up correctly. Much functionality unavailable."
                #return None
        #self.mol.Center()
        if not (axis == None or angle == None):
            self.rotate(axis,angle)
        if not vector == None:
            self.translate(vector)

    def duplicate(self, read_file=False):
        """
        Return a deep-copy of myself or re-read the original file.

        read_file: if True, the original file is read in again instead
                   of duplicating the molecule as it is.
        """
        if read_file:
            return read_from_file(self.fileinfo['name'],fileformat=self.fileinfo['format'],conf_nr=self.fileinfo['conf_nr'],ff=self.fileinfo['ff'])
        else:
            return molecule(self.mol,ff=self.ffname,fileinfo=self.fileinfo)

    
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
        print >>sys.stderr,"WARNING: the method 'optimize' heavily leaks memory and has to be rewritten to use only"
        print >>sys.stderr,"forcefields and no pybel. I'm not even sure whether it works properly."
        p_tempmol=p.Molecule(self.mol)
        p_tempmol.localopt(forcefield=self.ffname,steps=steps)

    def set_charge_method(self,method):
        self.charge_method = method

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
        vec = op.vector3(*axis)
        self.mol.Rotate(vec,angle)
        del vec
    
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
        vec = op.vector3(*vector)
        self.mol.Translate(vec)
        del vec

    def move_closer(self,part1,part2,stepsize=0.1,vdw_factor=0.9,vdw_added=0.0,vec=None):
        """
        Move two parts of an aggregate closer together. Indices start at 0.

        part1:      index indicating the first molecule in the aggregate that shall
                    be moved closer to another one
        part2:      second index
        sstepsize:  stepsize for movement (good value: 0.2)
        vdw_factor: factor by which all vdW-radii will be multiplied (default: 0.9)
        vdw_added:  value that is added to all vdW-radii (default: 0.0)
        """
        if vec is None:
            self.mol.MovePartsCloser(part1,part2,stepsize,vdw_factor,vdw_added)
        else:
            vtemp = _double_array(vec)
            self.mol.MovePartsCloser(vtemp, part1,part2,stepsize,vdw_factor,vdw_added)
            del vtemp
    
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
        vec = _double_array(vector)
        ax  = _double_array(axis)
        self.mol.AppendMolecule(mol.mol,vec,ax,angle)
        del vec,ax
    
    def write(self,filename,overwrite='False',fileformat='xyz'):
        """
        Write the data of the molecule to disk. 
    
        filename: name of the file INCLUDING the extension
        overwrite: shall the output file be overwritten or not
        fileformat: output file format (anything that openbabel can write)
        """
        if fileformat==None:
            fileformat=guess_format(filename)
        p.Molecule(self.mol).write(fileformat,filename,overwrite=overwrite)
    
    def align(self,point,main3,main2):
        """
        Align the last two main axes of a molecule to the two given axes and
        move the center to the given coordinate.
    
        point: 3-element sequence defining the new center of the molecule 
               (not mass weighed)
        main3: 3-element sequence defining the new 3rd main axis (longest extent)
        main2: 3-element sequence defining the new 2nd main axis
    
        """
        for vec in [[point,"point"],[main3,"third axis"],[main2,"second axis"]]:
            if not len(vec[0]) == 3:
                raise IndexError("Variable "+vec[1]+" not of the correct length, needs 3 elements not "+str(len(vec[0])))
        if (sum([abs(v) for v in main3])>0 and sum([abs(v) for v in main2])>0):
            poi = _double_array(point)
            ma3 = _double_array(main3)
            ma2 = _double_array(main2)
            self.mol.Align(poi,ma3,ma2)
            del poi,ma3,ma2

    def mirror(self,normal,point,center_it=False):
        """
        Mirror the molecule either by point inversion or by mirroring at a 
        plane.
        Align the first two main axes of a molecule to the two given axis and
        move the center to the given coordinate.
    
        normal: 3-element sequence defining the normal vector of the plane. If
                it is [0,0,0], point inversion will be performed.
        point:  3-element sequence defining either the inversion point or
                a point in the plane (Hessian normal form)
    
        """
        nor = _double_array(normal)
        poi = _double_array(point)
        self.mol.Mirror(nor,poi,center_it)
        del nor,poi
    
    def part_molecule(self,normal_vector,coordinate):
        """
        Get an openbabel molecule containing all those atoms that are one one side of a plane given
        in the Hessian normal form.
    
        normal_vector: normal vector of Hessian normal form (3-element list)
        coordinate: 3d-Cartesian coordinates of one point in the plane
        """
        tempmol=op.OBMol()
        nor = _double_array(normal_vector)
        coo = _double_array(coordinate)
        self.mol.PartMolecule(tempmol,nor,coo)
        del nor,coo
        return tempmol

    def part_molecule_mol(self,normal_vector,coordinate,side='left'):
        """
        Get a molecule containing all those atoms that are one one side of a plane given
        in the Hessian normal form.
    
        normal_vector: normal vector of Hessian normal form (3-element list)
        coordinate: 3d-Cartesian coordinates of one point in the plane
        side: If 'left', print those atoms on the side where the normal vector points.
              If anything else, print all those on the opposite side.
        """
        if side=='right':
            normal_vector=[-i for i in normal_vector]
        elif side=='left':
            pass
        else:
            raise ValueError("Side must be either left or right")
        return molecule(self.part_molecule(normal_vector,coordinate),ff=None)
    
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
        if side=='right':
            normal_vector=[-i for i in normal_vector]
        elif side=='left':
            pass
        else:
            raise ValueError("Side must be either left or right")
        p.Molecule(self.part_molecule(normal_vector,coordinate)).write(fileformat,filename,overwrite=overwrite)

    def get_partial_charges(self,method=None):
        """
        Return a list of all partial charges of the atoms
        according to the specified method. List has the
        same order as the atoms in the molecule.
        """
        #print >>sys.stderr,"WARNING: this function uses Gasteiger charges.\nWill be improved to also use other partitioning methods."
        if method is None:
            method = self.charge_method
        method = method.lower()
        tmp_charges = op.OBChargeModel.FindType(method)
        if tmp_charges is not None:
            if tmp_charges.ComputeCharges(self.mol):
                partialcharges = list(tmp_charges.GetPartialCharges())
            else:
                raise OpenBabelError("Error while partitioning partial charges.")
        else:
            raise ValueError("Method '"+method+"' is not a known method for partitioning partial charges. See 'obabel -L charges' for partitioning methods.")
        del tmp_charges
        #qeq does deliver charges of the opposite sign as the rest
        if method == 'qeq':
            partialcharges = [-q for q in partialcharges]
        return partialcharges

    def get_charges(self):
        """
        Return a list of all charges of the atoms according to ther element
        numbers. List has the same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        charges=[0.0]*self.mol.NumAtoms()
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            charges[idx-1] = a.GetAtomicNum()
        return charges

    def get_dipole_moment(self,method=None):
        """
        Return a list of all charges of the atoms according to ther element
        numbers. List has the same order as the atoms in the molecule.
        """
        charges=self.get_partial_charges(method=method)
        coordinates=self.get_coordinates()
        px=0
        py=0
        pz=0
        for (x,y,z),c in zip(coordinates,charges):
            px += x*c
            py += y*c
            pz += z*c
        return [px,py,pz]
        del charges, coordinates

    def get_coordinates(self):
        """
        Return a list of all cartesian coordinates of the atoms
        according to the current force field. List has the
        same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        coordinates=[None]*self.mol.NumAtoms()
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            coordinates[idx-1] = [a.GetX(),a.GetY(),a.GetZ()]
        return coordinates

    def get_center(self):
        """
        Return the non-mass-weighted center of the molecule
        """
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        return np.mean(np.array(self.get_coordinates()),axis=0)

    def get_main_axes(self):
        """
        Return the last 2 main axes of the molecule
        """
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        center = np.array(self.get_center())
        coords = np.array(self.get_coordinates())-center
        #mat is the tensor of inertia
        mat    = np.sum(np.array([ [[y*y+z*z,-x*y,-x*z],[-x*y,x*x+z*z,-y*z],[-x*z,-y*z,x*x+y*y]] for x,y,z in coords]),axis=0)
        eigvals,eigvecs = np.linalg.eig(mat)
        #the eigenvectors are stored in the coloumns of eigvecs
        #so it is transposed to have easy access to them
        eigvecs = -eigvecs.T
        main3,main2,main1 = sorted(zip(eigvals,eigvecs),key=lambda e: e[0])
        return main3[1],main2[1]

    def get_align_matrix(self,main3,main2):
        """
        Return the composite rotation matrix that would align the third and
        second main axes to the given axes.
        """
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        #c_ stands for current
        c_main3,c_main2 = self.get_main_axes()

        tempvec = np.cross(main3, c_main3)
        angle = _VectorAngle(main3, c_main3)
        mat1  = _RotMatrixAboutAxisByAngle(tempvec,angle)

        c_main2 = np.dot(mat1,c_main2)
        tempvec = np.cross(main2, c_main2)
        angle = _VectorAngle(main2, c_main2)
        mat2  = _RotMatrixAboutAxisByAngle(tempvec,angle)

        return np.dot(mat2,mat1)

    def get_vdw_radii(self):
        """
        Return a list of all van-der-Waals radii of the atoms
        according to the current force field. List has the
        same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        vdw_radii=[0.0]*self.mol.NumAtoms()
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            vdw_radii[idx-1] = op.etab.GetVdwRad(a.GetAtomicNum())
        return vdw_radii

    def get_names(self):
        """
        Return a list of all element symbolx radii of the atoms
        List has the same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        names=[]
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            names.append(op.etab.GetSymbol(a.GetAtomicNum()))
        return names

    def get_colours(self):
        """
        Return a list of all element symbolx radii of the atoms
        List has the same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        colours=[]
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            colours.append(op.etab.GetRGB(a.GetAtomicNum()))
        return colours

    def get_masses(self):
        """
        Return a list of all atomic masses of the atoms
        according to the current force field. List has the
        same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        masses=[0.0]*self.mol.NumAtoms()
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            masses[idx-1] = op.etab.GetMass(a.GetAtomicNum())
        return masses

    def get_vdw_surface_potential(self, nr_refinements=1, vertex='center', weights=None, triangulation=None, shrink_factor=0.95, charges=None, skip_potential=False):
        """
        Compute the static electric potential on a discretized van-der-Waals
        surface.

        nr_refinements: number of refinement steps for the skin surface.
                        The higher the number the more vertices it will have.
        vertex: are the vertices to be returned the centers of the triangles
                (default) or the corners (values 'center' and 'corners'). This
                is also where the potential will be computed.
        weights: (optional) give a list that will contain the surface area of
                 each triangle
        triangulation: (optional) a list that will contain the triangulation.
                       This is useful for 3d plotting
        shrink_factor: the shrink factor for the generation of the skin surface.
                       Must be >0 and <1. The bigger the tighter the surface will be.
        charges:       declare positions for localized charges and their values.
                       Format must be [[vec1,vec2,vec3],[c1,c2,c3]]
                       with ci being the charges and veci being 3-element vectors
        skip_potential:whether or not to skip the computation of the potential
                       If True, an emtpy list will be returned. This is useful
                       if you only want to have the triangulation

        Example in 2D for nr_points=12
        . : point on the sphere's surface
        X : center of the sphere

        Not overlapping => 12 points per sphere
            ...   ...
           .   . .   .
           . X . . X .
           .   . .   .
            ...   ...

        Overlapping => points that would be within the other sphere
                       are removed
                    => 9 points per "sphere"
            ......
           .      .
           . X  X .
           .      .
            ......
        """
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        if not supported["FireDeamon"][0]:
            raise MissingModuleError("Functionality requested that needs libFireDeamon but there was an error while importing the module.",supported["FireDeamon"][1])

        if charges == None:
            partialcharges=self.get_partial_charges()
            chargecoordinates=self.get_coordinates()
        else:
            chargecoordinates,partialcharges=charges

        vdw_radii=self.get_vdw_radii()
        coordinates=self.get_coordinates()
    
        lengths,face_indices,corners,normals = fd.SkinSurfacePy(shrink_factor,coordinates,vdw_radii,refinesteps=nr_refinements)
        triangles=[[np.array(corners[i]) for i in face] for face in face_indices]
        if weights!=None:
            trig_areas = [0.5*np.linalg.norm(np.cross(f[1]-f[0],f[2]-f[0])) for f in triangles]
            for p in trig_areas: weights.append(p)
        if triangulation!=None:
            for t in triangles: triangulation.append(t)

        if vertex=='center':
            trig_centres = [np.mean(f,axis=0) for f in triangles]
            if skip_potential:
                potential=[]
            else:
                potential = ep.potential_at_points(trig_centres, partialcharges, chargecoordinates)
            return trig_centres, potential
        elif vertex=='corners':
            if skip_potential:
                potential=[]
            else:
                potential = ep.potential_at_points(corners, partialcharges, chargecoordinates)
            return corners, potential
        else:
            raise ValueError("Wrong vertex type '"+vertices+"' specified.")

    def get_vdw_surface(self, get='faces', nr_refinements=1, shrink_factor=0.95, povray=0):
        """
        Conpute the discretized van-der-Waals surface.

        get: are the vertices to be returned the centers of the triangles
             the corner or shall the faces be returned (values 'center', 
             'corners' and 'faces'). Faces is the default.
        nr_refinements: number of refinement steps for the skin surface.
                        The higher the number the more vertices it will have.
        shrink_factor: the shrink factor for the generation of the skin surface.
                       Must be >0 and <1. The bigger the tighter the surface will be.
        povray: int, optional (default: 0)
            If >0, also to return the face indices and the bare vertex
            coordinates as a second and third list, respectively. This can be
            used to plot the surface using programmes such as PovRay. The
            resolution of the auto-generated PovRay plots will be the given
            value times the OpenGL resolution.

        Example in 2D for nr_points=12
        . : point on the sphere's surface
        X : center of the sphere

        Not overlapping => 12 points per sphere
            ...   ...
           .   . .   .
           . X . . X .
           .   . .   .
            ...   ...

        Overlapping => points that would be within the other sphere
                       are removed
                    => 9 points per "sphere"
            ......
           .      .
           . X  X .
           .      .
            ......
        """
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        if not supported["FireDeamon"][0]:
            raise MissingModuleError("Functionality requested that needs libFireDeamon but there was an error while importing the module.",supported["FireDeamon"][1])

        vdw_radii=self.get_vdw_radii()
        coordinates=self.get_coordinates()
    
        #lengths,face_indices,corners,normals = fd.SkinSurfacePy(shrink_factor,coordinates,vdw_radii,refinesteps=nr_refinements)
        lengths,face_indices,corners,normals = fd.SkinSurfacePy(shrink_factor,coordinates,vdw_radii,refinesteps=nr_refinements)
        triangles=[[np.array(corners[i]) for i in face] for face in face_indices]

        if get=='center':
            trig_centres = [np.mean(f,axis=0) for f in triangles]
            if povray>0:
                return trig_centres,face_indices,corners,normals
            else:
                return trig_centres
        elif get=='corners':
            if povray>0:
                return corners,face_indices,corners,normals
            else:
                return corners
        elif get=='faces':
            if povray>0:
                return triangles,face_indices,corners,normals
            else:
                return triangles 
        else:
            raise ValueError("Wrong vertex type '"+get+"' specified.")

    def get_bond_map(self,unique=True,no_hydrogen=False):
        """
        Produce a list of all bonds in a molecule as known by the current force field.

        unique: if True, give back an irreducible list of bonds in the form of tuples
                of indices
                if False, give back a complete list of bonds, i.e. every atom in a bond
                is once the first and once the second element in one of the tuples
        no_hydrogen: whether or not to exclude hydrogens from the list
        """
        bondmap=[]
        for bond_id in range(0,self.mol.NumBonds()):
            b=self.mol.GetBond(bond_id)
            if no_hydrogen and b.GetBeginAtom().IsHydrogen() or b.GetEndAtom().IsHydrogen():
                continue
            bondmap.append((b.GetBeginAtomIdx()-1,b.GetEndAtomIdx()-1))
        if unique:
            bondmap=[tuple(sorted(b,key=lambda x:x)) for b in bondmap]
            bondmap=sorted(list(set(bondmap)),key=lambda x:x[0]*(len(bondmap)+1)+x[1])
        else:
            bondmap=sorted(bondmap,key=lambda x:x[0]*(len(bondmap)+1)+x[1])
        return bondmap

    def visualize(self,zoom=1,align_me=True,point=[0.0,0.0,0.0],main3=[1,0,0],main2=[0,1,0],nr_refinements=1,method='simple',title="Molecule Visualization",resolution=(1024,768),high_contrast=False,spherescale=1,rendertrajectory=None,charges=None,potential=None,invert_potential=False,config=None,savefile=None,povray=0):
        """
        This function is a wrapper for visualizing the molecule using OpenGL.
        The molecule will be aligned prior to visualization.
        A helper module is being used. This has been done so that the lengthy
        visualization functions can reside in a different file.

        zoom: a zoom factor
        align_me: whether or not to align the molecule
                  prior to visualization
        point, main3, main2: see function "align"
        nr_refinements: number of subdivision steps after skin
                        surface generation
        method: if 'complex', visualize electrostatic potentials on
                a vdW-surface.
                If 'simple', visualize the atoms as coloured spheres.
        """
        try:
            from . import visualize_molecule as vm
        except ImportError as e:
            raise ImportError("Error importing helper module visualize_molecule",e)
        if method=='complex':
            if align_me:
                translate_before = -np.array(self.get_center())
                translate_after = np.array(point)
                rotate = self.get_align_matrix(main3,main2)
                manip_func = lambda e: np.dot(rotate,(np.array(e)+translate_before))+translate_after
            else:
                manip_func = None
            vm.PlotGL_Surface(self,zoom,nr_refinements=nr_refinements,title=title,resolution=resolution,high_contrast=high_contrast,rendertrajectory=rendertrajectory,charges=charges,ext_potential=potential,invert_potential=invert_potential,config=config,manip_func=manip_func,savefile=savefile,povray=povray)
        elif method=='simple':
            if align_me:
                self.align(point,main3,main2)
            vm.PlotGL_Spheres(self,zoom,title=title,resolution=resolution,spherescale=spherescale,rendertrajectory=rendertrajectory)
        else:
            raise ValueError("Selected method must be either complex or simple")

    def HLB_value(self, nr_refinements=0,anglestep=5,anglerange=10,no_hydrogen=True):
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        try:
            from . import hlb_value as hlb
        except ImportError as e:
            raise ImportError("Error importing helper module hlb_value",e)
        hlb_value,normal_vector,coordinate = hlb.get_HLB(self,nr_refinements,anglestep_degree=anglestep,anglerange_degree=anglerange,no_hydrogen=no_hydrogen)
        return hlb_value,normal_vector,coordinate

    def rmsd(self,molecule,print_result=False):
        if not supported["numpy"][0]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        selftemp=self.duplicate()
        othertemp=molecule.duplicate()
        selftemp.align([0,0,0],[1,0,0],[0,1,0])
        othertemp.align([0,0,0],[1,0,0],[0,1,0])
        selfcoords=np.array(selftemp.get_coordinates())
        othercoords=np.array(othertemp.get_coordinates())
        if not(selfcoords.shape == othercoords.shape):
            raise ValueError("The molecule to compare against does not have the same number of atoms.")
        diff=selfcoords-othercoords
        result_rmsd=np.sqrt(np.sum(diff*diff)/(3*len(selfcoords)))
        result_maxdeviation_single=np.max(np.abs(diff))
        result_maxdeviation_whole=np.max(np.linalg.norm(diff,axis=1))
        del selfcoords, othercoords
        del selftemp, othertemp
        result=(result_rmsd,result_maxdeviation_single,result_maxdeviation_whole)
        if print_result:
            print "RMSD: %.2e | Max deviation in one coordinate: %.2e | Max deviation for a single atom: %.2e"%result
        return result
