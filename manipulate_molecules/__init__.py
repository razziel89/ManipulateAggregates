"""
This script is able to arrange molecules arbitrarily.
Instead of doing so itself it rather provides the functions needed to do so.
Everything is wrapped in a class and uses openbabel data structures.
"""

global supported
supported={}

class ManipulateMoleculesError(Exception):
    pass

class NoOpenbabelError(ManipulateMoleculesError):
    pass

class MissingModuleError(ManipulateMoleculesError):
    pass

class WrongVertexError(ManipulateMoleculesError):
    pass

class WrongMethodError(ManipulateMoleculesError):
    pass

class WrongSideError(ManipulateMoleculesError):
    pass

class NotEnoughConformersError(ManipulateMoleculesError):
    pass

import re
try:
    import pybel as p
except ImportError as e:
    raise NoOpenbabelError("Pybel could not be imported. Please install openbabel with Python bindings.",e)

import sys
try:
    import numpy as np
    supported["numpy"]=(True,)
except ImportError as e:
    supported["numpy"]=(False,e)

try:
    import alphashapes as alpha
    supported["alpha"]=(True,)
except ImportError as e:
    supported["alpha"]=(False,e)
    
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

class FiletypeException(Exception):
    pass

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

def read_from_file(filename,fileformat=None,conf_nr=1):
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
    if conf_nr==1:
        return molecule(p.readfile(fileformat,filename).next().OBMol)
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
                raise NotEnoughConformersError("You requested conformer number %d but there are only %d present in the file."%(max(conf_nr),len(conformers)))
            return [molecule(conformers[i-1].OBMol) for i in conf_nr]
        else:
            if conf_nr>len(conformers):
                raise NotEnoughConformersError("You requested conformer number %d but there are only %d present in the file."%(conf_nr,len(conformers)))
            return molecule(conformers[conf_nr-1].OBMol)

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
#       op.OBGenericData.Clone(self.mol.OBGenericData,mol)
#       self.mol.CloneData(mol.GetData(op.OBGenericData_swigregister))
        if not ff in p.forcefields:
            return None
        self.ff=op.OBForceField.FindForceField(ff)
        self.ffname=ff
        if  self.ff.Setup(self.mol) == 0:
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
        if fileformat==None:
            fileformat=guess_format(filename)
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
        Get an openbabel molecule containing all those atoms that are one one side of a plane given
        in the Hessian normal form.
    
        normal_vector: normal vector of Hessian normal form (3-element list)
        coordinate: 3d-Cartesian coordinates of one point in the plane
        """
        tempmol=op.OBMol()
        self.mol.PartMolecule(tempmol,double_array(normal_vector),double_array(coordinate))
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
            raise WrongSideError("Side must be either left or right")
        return molecule(self.part_molecule(normal_vector,coordinate))
    
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
            raise WrongSideError("Side must be either left or right")
        p.Molecule(self.part_molecule(normal_vector,coordinate)).write(fileformat,filename,overwrite=overwrite)

    def get_partial_charges(self):
        """
        Return a list of all partial charges of the atoms
        according to the current force field. List has the
        same order as the atoms in the molecule.
        """
        a=op.OBAtom()
        partialcharges=[0.0]*self.mol.NumAtoms()
        for idx in range(1,self.mol.NumAtoms()+1):
            a = self.mol.GetAtom(idx)
            partialcharges[idx-1] = a.GetPartialCharge()
        return partialcharges

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

    def get_vdw_surface_potential(self, nr_refinements=1, vertex='center', weights=None, triangulation=None, shrink_factor=0.95):
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

        partialcharges=self.get_partial_charges()
        coordinates=self.get_coordinates()
        vdw_radii=self.get_vdw_radii()
    
        lengths,face_indices,corners = fd.SkinSurfacePy(shrink_factor,coordinates,vdw_radii,refinesteps=nr_refinements)
        triangles=[[np.array(corners[i]) for i in face] for face in face_indices]
        if weights!=None:
            trig_areas = [0.5*np.linalg.norm(np.cross(f[1]-f[0],f[2]-f[0])) for f in triangles]
            for p in trig_areas: weights.append(p)
        if triangulation!=None:
            for t in triangles: triangulation.append(t)

        if vertex=='center':
            trig_centres = [np.mean(f,axis=0) for f in triangles]
            potential = ep.potential_at_points(trig_centres, partialcharges, coordinates)
            return trig_centres, potential
        elif vertex=='corners':
            potential = ep.potential_at_points(corners, partialcharges, coordinates)
            return corners, potential
        else:
            raise WrongVertexError("Wrong vertex type '"+vertices+"' specified.")

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

    def visualize(self,zoom=1,align_me=True,point=[0.0,0.0,0.0],main1=[0,0,1],main2=[0,1,0],nr_refinements=1,method='complex',title="Molecule Visualization",resolution=(1024,768),high_contrast=False):
        """
        This function is a wrapper for visualizing the molecule using OpenGL.
        The molecule will be aligned prior to visualization.
        A helper module is being used. This has been done so that the lengthy
        visualization functions can reside in a different file.

        zoom: a zoom factor
        align_me: whether or not to align the molecule
                  prior to visualization
        point, main1, main2: see function "align"
        nr_refinements: number of subdivision steps after skin
                        surface generation
        method: if 'complex', visualize electrostatic potentials on
                a vdW-surface.
                If 'simple', visualize the atoms as coloured spheres.
        """
        if align_me:
            self.align(point,main1,main2)
        try:
            from . import visualize_molecule as vm
        except ImportError as e:
            raise ImportError("Error importing helper module visualize_molecule",e)
        if method=='complex':
            vm.PlotGL_Surface(self,zoom,nr_refinements=nr_refinements,title=title,resolution=resolution,high_contrast=high_contrast)
        elif method=='simple':
            vm.PlotGL_Spheres(self,zoom,title=title,resolution=resolution)
        else:
            raise WrongMethodError("Selected method must be either complex or simple")

    def HLB_value(self, nr_refinements=0,anglestep=5,anglerange=10,no_hydrogen=True):
        if not supported["numpy"]:
            raise MissingModuleError("Functionality requested that needs numpy but there was an error while importing the module.",supported["numpy"][1])
        try:
            from . import hlb_value as hlb
        except ImportError as e:
            raise ImportError("Error importing helper module hlb_value",e)
        hlb_value,normal_vector,coordinate = hlb.get_HLB(self,nr_refinements,anglestep_degree=anglestep,anglerange_degree=anglerange,no_hydrogen=no_hydrogen)
        return hlb_value,normal_vector,coordinate