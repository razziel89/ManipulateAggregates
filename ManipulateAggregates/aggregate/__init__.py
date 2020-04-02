"""Definitions of the agg class and some auxilliary function.

If an external module cannot be imported, all the functionality that does not
require this module is still supported.
"""

# This file is part of ManipulateAggregates.
#
# Copyright (C) 2016 by Torsten Sachse
#
# ManipulateAggregates is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ManipulateAggregates is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import re
import sys
import copy
import itertools
import logging

logger = logging.getLogger(__name__)

import maagbel
from maagbel import pybel

from . import hlb
from . import visualize

from .. import orbitalcharacter

global SUPPORTED
# remember which submodules could successfully be imported and which not
SUPPORTED = {}


class ManipulateMoleculesError(Exception):
    """Base error class"""

    pass


E_UNIT_CONVERSION = {
    "kJ/mol->meV": 0.09648500,
    "kcal/mol->meV": 0.02306035,
    "meV->kJ/mol": 1.0 / 0.09648500,
    "meV->kcal/mol": 1.0 / 0.02306035,
    "kJ/mol->kcal/mol": 0.09648500 / 0.02306035,
    "kcal/mol->kJ/mol": 0.02306035 / 0.09648500,
}

## default meshing criteria for the iso surface generation
global MESH_CRITERIA
MESH_CRITERIA = [5.0, 0.2, 0.2]


class MissingModuleError(ManipulateMoleculesError):
    """Raised if a required module could not be imported."""

    pass


class FiletypeException(ManipulateMoleculesError):
    """Raised if an auto-determined file type is unknown."""

    pass


class OpenBabelError(ManipulateMoleculesError):
    """Raised if a computation in OpenBabel failed."""

    pass


try:
    import numpy

    SUPPORTED["numpy"] = (True,)
except ImportError as e:
    SUPPORTED["numpy"] = (False, e)

try:
    import FireDeamon as fd

    SUPPORTED["FireDeamon"] = (True,)
except ImportError as e:
    SUPPORTED["FireDeamon"] = (False, e)

try:
    from ..collection import read as fdread

    SUPPORTED["ManipulateAggregates.collection.read"] = (True,)
except ImportError as e:
    SUPPORTED["ManipulateAggregates.collection.read"] = (False, e)


def _assert_supported(key):
    """Make sure that a module that is required is actually supported on this machine.

    Args:
        key: (string) module name

    Raises:
        ManipulateAggregates.aggregate.MissingModuleError.
    """
    value = SUPPORTED.get(key, None)
    if value is not None:
        if not value[0]:
            raise MissingModuleError(
                "Functionality requested that needs %s but there was an error while importing the module:"
                % (key),
                value[1],
            )
    else:
        raise ValueError("Module %s unknown to the import checker." % (key))


global FILETYPEDICT
# copy all formats known to maagbel into a new dictionary
FILETYPEDICT = {}
for entry in pybel.informats:
    FILETYPEDICT[entry] = entry
# FILETYPEDICT={entry:entry for entry in pybel.informats}
# add some custom entries
FILETYPEDICT["mop"] = "mopin"


def _double_array(mylist):
    """Create a C array of doubles from a list."""
    c = maagbel.doubleArray(len(mylist))
    for i, v in enumerate(mylist):
        c[i] = v
    return c


def _vector3(mylist):
    if len(mylist) != 3:
        raise ValueError("vector3 has to have exactly 3 elements.")
    return maagbel.vector3(*mylist)


def _obbitvec(mask):
    """Get a maagbel bit vector representation of the mask.

    Returns:
        a swig proxy to the maagbel bitvector mask
    """
    if True in (i < 0 for i in mask):
        raise ValueError("Have to start with the zeroth entry.")
    bvmask = maagbel.OBBitVec()
    for i in mask:
        bvmask.SetBitOn(i)
    return bvmask


def guess_format(filename, no_except=False):
    """Try to guess the file format.
    
    If the filename contains no dot, try to match the whole filename against
    OpenBabel's filetype database. Useful for e.g. CONTCAR files that seldomnly
    contain a dot.

    Args:
        filename: (string) the name of the file whose type to guess
        no_except: (bool) if False, an exception is raised if the filetype is not known.
            If True, only None will be returned in that case.

    Returns:
        the file type or None (depending on no_except)

    Raises:
        ManipulateAggregates.aggregate.FiletypeException.
    """
    temp = filename.split(".")
    if len(temp) == 1:
        extension = filename
    else:
        extension = temp[-1]
    filetype = FILETYPEDICT.get(extension, None)
    if filetype is None and not no_except:
        raise FiletypeException(
            "Filetype of file " + filename + " not known to maagbel.", e
        )
    return filetype


def read_from_file(filename, fileformat=None, conf_nr=1, ff="mmff94"):
    """Convert content of file to ManipulateAggregates.aggregate.agg.

    Guesses filetype if none specified. Also converts "~" to the proper home
    directory.

    Args:
        filename: (string) path to the file
        fileformat: (string or None) guess type if None, otherwise use specified
            filetype
        conf_nr: (int or iterable of ints) index/indices of those conformers to be
            loaded. Special keyword 'all' will return all conformers in file.  Special
            keywords "first" and "last" return the first and last of the conformers,
            respectively
        ff: (string) name of the force field to use. Run "manipagg --list forcefields"
            to get supported ones.

    Returns:
        object of ManipulateAggregates.aggregate.agg

    Raises:
        ValueError.
    """
    if fileformat is None:
        fileformat = guess_format(filename)
    if re.match("~.*" + os.sep, filename):
        homedir = filename.split(os.sep)[0]
        if re.match("^~$", homedir):
            filename = re.sub("^~", os.environ["HOME"], filename)
        else:
            username = homedir.split("~")[1]
            homedir = os.sep.join(os.environ["HOME"].split(os.sep)[:-1])
            filename = re.sub("^~", homedir + os.sep, filename)

    info = {"name": filename, "format": fileformat, "conf_nr": conf_nr, "ff": ff}

    if conf_nr not in ("all", "first", "last") and conf_nr < 1:
        raise ValueError("Conformers are counted starting at 1.")
    elif conf_nr == "last":
        obagg = None
        for m in pybel.readfile(fileformat, filename):
            obagg = m
        # obagg = reduce(lambda x,y: y, (m for m in pybel.readfile(fileformat,filename))).OBMol
        return agg(obagg.OBMol, ff=ff, info=info)
    else:
        # conf_nr_iter = itertools.count(start=0,step=1)
        conf_nr_iter = itertools.count()
        if conf_nr == "first":
            conf_nr = 1
            conf_nr_match = lambda n: n == 0
            max_conf_nr = conf_nr
        elif isinstance(conf_nr, int):
            conf_nr_match = lambda n: n == conf_nr - 1
            max_conf_nr = conf_nr
        elif conf_nr == "all":
            conf_nr_match = lambda n: True
            max_conf_nr = float("inf")
        else:
            raise ValueError(
                "Wrong type for argument conf_nr, type is: %s" % (str(type(conf_nr)))
            )

        def gen_conformers():
            for pmol in pybel.readfile(fileformat, filename):
                conf_count = next(conf_nr_iter)
                if conf_count >= max_conf_nr:
                    return
                else:
                    yield (conf_count, pmol.OBMol)

        conformers = list(gen_conformers())
        if len(conformers) > 0:
            first = True
            for ccount, cagg in conformers:
                if conf_nr_match(ccount):
                    if first:
                        obagg = cagg
                    else:
                        obagg.AddConformer(cagg.GetCoordinates(), True)
            return agg(obagg, ff=ff, info=info)
        else:
            if isinstance(conf_nr, int):
                raise ValueError(
                    "You requested conformer number %d but there are fewer in the file."
                    % (conf_nr,)
                )
            else:
                return None


def _RotMatrixAboutAxisByAngle(axis, angle):
    """
    Taken from maagbel from file matrix3x3.cpp, method RotAboutAxisByAngle.
    Generate a rotation matrix about an arbitrary axis by an arbitrary angle.
    Angle has to be in radians.
    """
    mat = numpy.identity(3, dtype=float)
    theta = angle
    s = numpy.sin(theta)
    c = numpy.cos(theta)
    t = 1.0 - c

    vtmp = numpy.array(axis, dtype=float)
    if not (len(vtmp.shape) == 1 and vtmp.shape[0] == 3):
        raise ValueError(
            "Given axis must have shape (3,) but it has shape " + str(vtmp.shape)
        )
    if numpy.linalg.norm(vtmp) > 0.001:
        vtmp /= numpy.linalg.norm(vtmp)

        x, y, z = vtmp

        mat[0][0] = t * x * x + c
        mat[0][1] = t * x * y + s * z
        mat[0][2] = t * x * z - s * y

        mat[1][0] = t * y * x - s * z
        mat[1][1] = t * y * y + c
        mat[1][2] = t * y * z + s * x

        mat[2][0] = t * z * x + s * y
        mat[2][1] = t * z * y - s * x
        mat[2][2] = t * z * z + c

    return mat


def _VectorAngle(vector1, vector2):
    """
    Return the angle between two vectors in radians.
    """
    v1 = numpy.array(vector1)
    v2 = numpy.array(vector2)
    v1 /= numpy.linalg.norm(v1)
    v2 /= numpy.linalg.norm(v2)
    dp = numpy.dot(v1, v2)
    if dp < -1.0:
        dp = -1.0
    elif dp > 1.0:
        dp = 1.0
    return numpy.arccos(dp)


class dummy_ff:
    """Dummy class to catch cases when OpenBabel's force field could not be set up."""

    def Setup(self, mol):
        return False

    def GetUnit(self, mol):
        raise RuntimeError("Not to be called, only a dummy.")

    def SteepestDescent(self, mol):
        raise RuntimeError("Not to be called, only a dummy.")

    def GetCoordinates(self, mol):
        raise RuntimeError("Not to be called, only a dummy.")

    def Energy(self, mol):
        raise RuntimeError("Not to be called, only a dummy.")


class agg:
    """The agg class (agg stands for aggregate)

    This class is:
      1. a Python wrapper for the OBAggregate class of the modified OpenBabel
         version providing pythonic methods that are not very easy to use via
         OpenBabel's language bindings
      2. an extension to the above in the sense that, by interfacing with the other
         submodules of ManipulateAggregates, additional functionality is provided,
         such as:

            - estimation of the HLB value
            - visualization of an electrostatic potential on a molecular surface

      3. the base class used for the aggregate geometry estimation implemented in
         the submodule ManipulateAggregates.energyscan

    Attributes:
        obmol: (of type OBAggregate) the underlying OpenBabel data structure
        info: (dictionary) contains information about the file this object was created
            from
        vs: (dictionary) contains information about desired visualizations and surface
            generation
        cp: (dictionary) contains information about how to obtain this molecule's
            charges and potentials
        __internal__: (dictionary) internal data, not for direct use
        ff: (OBForceField) the underlying OpenBabel force field data structure
    """

    # Default config options for surface generation and visualization.
    #
    # See comments for explanations of each parameter.
    default_vs = {
        # The type of plot. "vdw" (van-der-Waals surface) and "iso" (isosurface) will
        # show the elctrostatic potential on these surfaces. The type "simple" shows
        # coloured vdW spheres
        "type": "simple",  # "simple", "vdw" or "iso"
        # The isovalue for type "iso". For electron densities, the default is a good
        # value.
        "isovalue": 0.005,  # float
        # The type of file from which to get the electron density.
        "isofiletype": "dx",  # only "dx" supported
        # The path to the file from which to get the electron density.
        "isofile": "",  # string
        # Criteria for generation of the iso surfaces. For further info, see docstring
        # of funtion IsosurfacePy in external module FireDeamon for further explanations
        # or the link http://doc.cgal.org/latest/Surface_mesher/index.html. The default
        # should be fine.
        "mesh_criteria": MESH_CRITERIA,  # list of 3 floats
        # Precision value used to compute the isosurface. A lower value results in more
        # highly discretized surfaces.
        "rel_precision": 1.0e-06,  # float >0
        # Which atoms shall be used for the generation of the iso surface. Please see
        # the keyword "atoms" for the method get_iso_surface for more information.
        "iso_atoms": "auto",  # int, list of ints or "all" or "noH" or "auto"
        # The number of refinement steps for the generation of vdW surfaces. A higher
        # value creates a more highly discretized surface but needs overproportionally
        # more memory. More than 3 should never be required. See get_vdw_surface.
        "refine": 1,  # int >0
        # The shrink factor during surface generation. The default corresponds to vdW
        # surfaces. See get_vdw_surface for more information.
        "shrink_factor": 0.95,  # float >0 and <1
        # Scale all vdW radii by this value prior to visualization. Used for
        # visualization types "simple" and "vdw".
        "vdw_scale": 1.0,  # float >0
        # How the colour scaling shall be obtained. The value "independent" causes
        # positive and negative colour scales to be independent. The value "dependent"
        # causes the overall maximum (in terms of the absolute value) to correspond to
        # the highest potential. If a tuple of 2 strings is given, the first is taken as
        # a regular expression and the second one as a directory. The color scale will
        # be chosen in such a way that the color scales defined in all the saved
        # visualization state files in the given directory (whose names match the
        # regular expression) are contained within the new color scale.
        "colorscale": "independent",  # "independent", "dependent" or "REGEX","DIR"
        # Scale the position of all surface vertices by this value.
        "zoom": 1.0,  # float >0
        # The title of the created OpenGL window.
        "title": "Molecule Visualization",  # string
        # The resolution of the OpenGL window in pixels. Also the base resolution of
        # images rendered by PoVRay.
        "resolution": (1024, 768),  # tuple of 2 ints >0
        # If True, red, black and blue correspond to a positive, vanishing and negative
        # electrostatic potential. If False, another color scheme is used (deep blue,
        # light blue, green, yellow and red) and deep blue/red correspond to a highly
        # positive/negative and green to a vanishing potential.
        "high_contrast": True,  # bool
        # A trajectory tells the script how to manipulate the visualization of the given
        # aggregate.  All images that are rendered can be saved to disk. The format is
        # comples. Please see the --render-help message of the "manipagg" executable
        # contained in this bundle.
        "renderpath": None,  # string
        # Whether or not the OpenGL window shall be shown. Only makes sense with a
        # renderpath.
        "hide": False,  # bool
        # If <=0, PoVRay support is switched off. If >0, the resolution of the OpenGL
        # plot will be scaled by this value to get the resolution for the corresponding
        # PoVRay image.
        "povray": 1,  # int >0
        # A file name (no path!!!) to which the visualization state can be saved. No
        # default.  Information about when and why this state was saved might be
        # prefixed to the name. Saving the visualization state only works with types
        # "vdw" and "iso".
        "savefile": None,  # string
        # Whether or not to save the visualization state at the beginning of the
        # visualization.
        "savestart": False,  # bool
        # Whether or not to save the visualization state at the end of the
        # visualization.
        "saveend": False,  # bool
        # Wether or not to align the aggregate prior to visualization. This is
        # recommened to make the aggregate better visible. The actual object will not be
        # changed, only the visualization will be adjusted.
        "align": True,  # bool
        # The non-mass-weighted center of the aggregate will be pot here.
        "align_center": [0.0, 0.0, 0.0],  # list of 3 floats
        # The third main axis (usually the longest extend) of the aggregate will point
        # in this direction.
        "align_main3": [1.0, 0.0, 0.0],  # list of 3 floats
        # The second main axis of the aggregate will point in this direction.
        "align_main2": [0.0, 1.0, 0.0],  # list of 3 floats
        # A filename to which a graphical representation of the color scale might be
        # saved (SVG file)
        "svgscale": None,
        # A rotation matrix that will be applied to all normal vectors prior to PoVRay
        # visualization
        "visrotmat": ([-1, 1, 0], 45),
    }
    # Default config options for obtaining charges and electrostatic potentials
    #
    # See comments for explanations of each parameter.
    default_cp = {
        # How the electrostatic potential shall be obtained. The value "empirical" uses
        # empirical methods implemented in OpenBabel (generation of partial charges,
        # which are then used to compute the potential) (fastest method). The value
        # "orbitals" uses information from quantum mechanical orbitals to compute the
        # potential (most accurate method). The value "interpolation" interpolates a
        # given electrostatic potential onto new points in space. The value "charge"
        # uses externally provided charges instead of, e.g., empirically determined ones
        # and then computes the potential.
        "type": "empirical",  # "empirical", "orbitals", "interpolation" or "charges"
        # How empirical charges shall be obtained. Everything supported by OpenBabel is
        # supported.
        "method": "mmff94",  # string
        # Whether or not charges obtained using empirical methods are partial charges or
        # not.
        "partial": True,  # bool
        # The type of file in which the quantum mechanical orbital data is stored.
        "orbfiletype": "molden",  # only "molden" supported so far
        # The path to the file in which the quantum mechanical orbital data is stored.
        "orbfile": "",  # string
        # The type of file in which the charge data is stored.
        "chargefiletype": "xyz",  # can be "xyz", "dx" or "cube"
        # The path to the file in which the charge data is stored.
        "chargefile": "",  # string
        # The type of file in which the potential data is stored.
        "potfiletype": "dx",  # can be "dx" or "xyz"
        # The path to the file in which the potential data is stored.
        "potfile": "",  # string
        # Whether or not the computed electrostatic potential should be inverted.
        "invert_potential": False,  # bool
        # How external ptoential data shall be interpolation. You can use "nearest" for
        # nearest neighbour interpolation or "distance" for inverse distance
        # interpolation. The latter uses the formula weight=1/(R-Norm())^E. See below.
        "interpolation": "nearest",  # can be "distance" or "nearest"
        # This is "E" in the above formula. Anything below 3 causes far-away values to
        # dominate.
        "int_exponent": 3,  # int >0
        # This is "R" in the above formula. The value 2 is the Eukledian norm.
        "int_root": 2,  # int >0
        # If the distance between a new and an old point is farther than this, consider
        # its weight to be zero. When computing the density from orbitals, also use this
        # cutoff for the distance between an atomic center and the point.
        "cutoff": 7.0,  # float >0
        # The total charge of this molecule. Only considered for "potfiletype"=="dx",
        # "chargefiletype"=="dx" and "chargefiletype"=="cube".
        "total_charge": 0,  # float
        # The property to compute. Used externally only.
        "property": "potential",  # can be "potential" or "density"
    }

    def __init__(self, obmol, ff="mmff94", info={}):
        """Constructor.

        Args:
            obmol: (OBAggregate) the underlying OpenBabel data structure, will be copied
            ff: (string) declare the forcefield associated with the molecule.  Needed to
                get the energy and perform a simple forcefield geometry optimization.
                Run "manipagg --list forcefield" to get supported ones. Can also be None
                (switched off)
            info: (dictionary) has keys 'name' and 'format' detailing the filename and
                format, respectively. Also, the keys 'conf_nr' and 'ff' for the used
                conformer number and force field are required.
        """
        self.obmol = maagbel.OBAggregate(obmol)
        self.info = copy.deepcopy(info)
        self.info["outformat"] = self.info.get("format", "nul")
        self.cp = copy.deepcopy(agg.default_cp)
        self.vs = copy.deepcopy(agg.default_vs)
        self.__internal__ = {
            "orb": None,
            "orbcfg": None,
            "pot": None,
            "potcfg": None,
            "cha": None,
            "chacfg": None,
        }
        if ff is not None:
            if not ff in pybel.forcefields:
                print("Force field not known to maagbel.", file=sys.stderr)
                self.ff = dummy_ff()
                self.info["ff"] = None
            else:
                self.ff = maagbel.OBForceField.FindForceField(ff)
                if self.ff.Setup(self.obmol) == 0:
                    print(
                        "Force field could not be set-up correctly. Much functionality unavailable.",
                        file=sys.stderr,
                    )
                    self.info["ff"] = None
                    self.ff = dummy_ff()
                else:
                    self.info["ff"] = ff
        else:
            self.info["ff"] = None
            self.ff = dummy_ff()

    def __del__(self):
        """Destructor."""
        del self.obmol
        del self.info
        del self.cp
        del self.vs
        del self.__internal__
        if self.ff is not None:
            del self.ff

    def duplicate(self, read_file=False):
        """Duplicate myself or re-read the original file.

        Please be aware that the constructor is called again. The __internal__
        dictionary is only shallowly copied, all the others are deep copies.

        Args:
            read_file: (bool) if True, the original file is read in
                again instead of duplicating the molecule as it is.

        Returns:
            object of ManipulateAggregates.aggregate.agg
        """
        if read_file:
            return read_from_file(
                self.info["name"],
                fileformat=self.info["format"],
                conf_nr=self.info["conf_nr"],
                ff=self.info["ff"],
            )
        else:
            newagg = agg(self.obmol, ff=self.info["ff"], info=self.info)
            newagg.cp = copy.deepcopy(self.cp)
            newagg.vs = copy.deepcopy(self.vs)
            newagg.__internal__ = copy.copy(self.__internal__)
            return newagg

    def get_pointgroup(self):
        """Determine the point group of the current aggregate.

        Returns:
            the name of the pointgroup.
        """
        sym = maagbel.OBPointGroup()
        sym.Setup(self.obmol)
        pg = sym.IdentifyPointGroup(0.01)
        del sym
        return pg

    def get_energy(self, unit="meV"):
        """Get the energy associated with the current geometry for the current forcefield.

        Units can be specified. Supported units are "kJ/mol", "kcal/mol" and "meV".

        Returns:
            energy in the specified units.
            
        """
        success = self.ff.Setup(self.obmol)
        if not success:
            raise OpenBabelError("Error setting up forcefield.")
        ffunit = self.ff.GetUnit()
        if ffunit != unit:
            try:
                conversion = E_UNIT_CONVERSION[ffunit + "->" + unit]
            except KeyError as e:
                raise ValueError(
                    "Unknown target (%s) or origin (%s) unit. I know: 'kJ/mol', 'kcal/mol' and 'meV'."
                    % (unit, ffunit),
                    e,
                )
        else:
            conversion = 1.0

        return self.ff.Energy(False) * conversion

    def optimize(self, steps=500):
        """Perform a sinple geometry optimization using the current forcefield.

        Args:
            steps: (int) number of optimization steps
        """
        success = self.ff.Setup(self.obmol)
        if not success:
            raise OpenBabelError("Error setting up forcefield.")
        self.ff.SteepestDescent(steps)
        self.ff.GetCoordinates(self.obmol)

    def set_cp(self, key, value):
        """Set an arbitrary configuration option of the dictionary for obtaining charges and potentials.

        If key and value are lists , key-value pairs will be assigned.
        If the lists are of unequal lenghts (e.g., n and n+m), only the
        first n key-value pairs will be treated.

        Args:
            key: (string) the property (or properties) to set
            value: (appropriate) the data associated with the key(s). See
                ManipulateAggregates.aggregate.agg.default_cp for possible
                options and values.
        """
        if isinstance(key, list) and isinstance(value, list):
            kviter = iter(zip(key, value))
            iterable = True
        else:
            iterable = False
        if iterable:
            for k, v in kviter:
                if k in agg.default_cp:
                    self.cp[k] = v
                else:
                    print("Key '%s' unknown, skipping." % (k), file=sys.stderr)
        else:
            if key in agg.default_cp:
                self.cp[key] = value
            else:
                print("Key '%s' unknown, skipping." % (key), file=sys.stderr)

    def set_vs(self, key, value):
        """Set an arbitrary configuration option of the dictionary for visualizations and surface generation.

        If key and value are lists, key-value pairs will be assigned.
        If the lists are of unequal lenghts (e.g., n and n+m), only the
        first n key-value pairs will be treated.

        Args:
            key: (string or list of strings) the property (or properties) to set
            value: (appropriate) the data associated with the key(s). See
                ManipulateAggregates.aggregate.agg.default_vs for possible
                options and values.
        """
        if isinstance(key, list) and isinstance(value, list):
            kviter = iter(zip(key, value))
            iterable = True
        else:
            iterable = False
        if iterable:
            for k, v in kviter:
                if k in agg.default_vs:
                    self.vs[k] = v
                else:
                    print("Key '%s' unknown, skipping." % (k), file=sys.stderr)
        else:
            if key in agg.default_vs:
                self.vs[key] = value
            else:
                print("Key '%s' unknown, skipping." % (key), file=sys.stderr)

    def get_cp(self, key):
        """Get an arbitrary configuration option of the dictionary that configures how charges and potentials are obtained.

        Args:
            key: (string) the property whose value you want to get

        Returns:
            the value associated with key or None if the key is not present
        """
        return self.cp.get(key, None)

    def get_vs(self, key):
        """Get an arbitrary configuration option of the dictionary for visualizations and surface generation.

        Args:
            key: (string) the property whose value you want to get

        Returns:
            the value associated with key or None if the key is not present
        """
        return self.vs.get(key, None)

    def set_bondlength(self, idx1, idx2, length, fix=None):
        """Adjust the length of a bond.
        
        If the bond connects two parts of a molecule that are otherwise not
        connected, those parts are moved with the respective atom. Otherwise,
        move only the 2 given atoms. There does not actually have to be a bond
        between the given atoms.
    
        Example original geometry

        >>>  _   _
        >>> |_|-|_|

        Adjust the middle bond to 3 times it's original length

        >>>  _     _
        >>> |_|---|_|

        Please note how both squares were moved along with the atoms
        comprising the bond.
    
        Args:
            idx1: (int) number of first atom that defines the bond
            idx2: (int) number of second atom that defines the bond
            length: (float) the new bond length (in Angstroms)
            fix: (1 or 2 or None) keep the first or second atom fixed
                and move only the other. if it is None, move both by half the required
                distance
        """
        bond = maagbel.OBBond()
        bond = self.obmol.GetBond(int(idx1), int(idx2))
        if fix is None:
            bond.SetLength(float(length))
        elif fix == 1:
            bond.SetLength(self.obmol.GetAtom(int(idx1)), float(length))
        elif fix == 2:
            bond.SetLength(self.obmol.GetAtom(int(idx2)), float(length))
        else:
            raise ValueError("The kwarg 'fix' has to be in (1,2,None)")

    def get_bondlength(self, idx1, idx2, projection=None):
        """Get the length of a bond.
        
        There does not actually have to be a bond between the given atoms. If
        projection is not None, the bond will be projected onto the given
        vector and the length of that projection will be given.
    
        Args:
            idx1: (int) number of first atom that defines the bond
            idx2: (int) number of second atom that defines the bond
            projection: (list of 3 floats) projection vector

        Returns:
            (possibly projected) bond length
        """
        a1 = self.obmol.GetAtom(int(idxa))
        a2 = self.obmol.GetAtom(int(idx2))
        pos1 = [a1.GetX(), a1.GetY(), a1.GetZ()]
        pos2 = [a2.GetX(), a2.GetY(), a2.GetZ()]
        if projection is None:
            return (
                (a1.GetX() - a2.GetX()) ** 2
                + (a1.GetY() - a2.GetY()) ** 2
                + (a1.GetZ() - a2.GetZ()) ** 2
            ) ** 0.5
        else:
            if len(projection) != 3:
                raise ValueError("Length of projection vector unequal 3.")
            pnorm = pow(sum(1.0 * p * p for p in projection), 0.5)
            return (
                (projection[0] * (a1.GetX() - a2.GetX()) / pnorm) ** 2
                + (projection[1] * (a1.GetY() - a2.GetY()) / pnorm) ** 2
                + (projection[2] * (a1.GetZ() - a2.GetZ()) / pnorm) ** 2
            ) ** 0.5

    def set_angle(self, idx1, idx2, idx3, angle):
        """Set the bond angle in deg.
        
        If the angle connects two parts of a molecule that are otherwise not
        connected, those parts are moved with the respective atom. See
        ManipulateAggregates.aggregate.agg.set_bondlength for a graphical example.
    
        Args:
            idx1: (int) number of first atom that defines the angle
            idx2: (int) number of second atom that defines the angle
            idx3: (int) number of third atom that defines the angle
            angle: (float) the new angle (in degrees)
        """
        self.obmol.SetAngle(int(idx1), int(idx2), int(idx3), float(angle))

    def get_angle(self, idx1, idx2, idx3):
        """Get the bond angle in deg. 
    
        Args:
            idx1: (int) number of first atom that defines the angle
            idx2: (int) number of second atom that defines the angle
            idx3: (int) number of third atom that defines the angle

        Returns:
            the angle defined by the three atoms in degrees
        """
        return self.obmol.GetAngle(int(idx1), int(idx2), int(idx3))

    def set_dihedral(self, idx1, idx2, idx3, idx4, angle):
        """Set the dihedral angle in deg.

        If the angle connects two parts of a molecule that are otherwise not
        connected, those parts are moved with the respective atom. See
        ManipulateAggregates.aggregate.agg.set_bondlength for a graphical example.
    
        Args:
            idx1: (int) number of first atom that defines the dihedral angle
            idx2: (int) number of second atom that defines the dihedral angle
            idx3: (int) number of third atom that defines the dihedral angle
            idx4: (int) number of fourth atom that defines the dihedral angle
            angle: (float) the new angle (in degrees)
        """
        self.obmol.SetDihedralAngle(
            int(idx1), int(idx2), int(idx3), int(idx4), float(angle)
        )

    def get_dihedral(self, idx1, idx2, idx3, idx4):
        """Get the dihedral angle in deg. 
    
        Args:
            idx1: (int) number of first atom that defines the dihedral angle
            idx2: (int) number of second atom that defines the dihedral angle
            idx3: (int) number of third atom that defines the dihedral angle
            idx4: (int) number of fourth atom that defines the dihedral angle

        Returns:
            the dihedral angle defined by the four atoms in degrees
        """
        return self.obmol.GetTorsion(int(idx1), int(idx2), int(idx3), int(idx4))

    def rotate(self, axis, angle, part=None):
        """Rotate the molecule around an axis by an angle.
    
        Args:
            axis: (list of 3 floats) rotate the geometry around this axis
            angle: (float) the angle for the rotation (in degrees)
            part: (None or int>=0) if an integer is provided, only treat
                that part. Please see ManipulateAggregates.aggregate.agg.tag_parts for a
                definition of "part".
        """
        vec = _double_array(axis)
        if part is None:
            self.obmol.Rotate(vec, float(angle))
        else:
            self.obmol.RotatePart(int(part), vec, float(angle))
        del vec

    def rotate_main(self, axis_index, angle, part=None):
        """Rotate the molecule around one of its main axis by an angle.
    
        Args: 
            axis_index: (int, 1, 2 or 3) the index of the main axis to rotate around
            angle: (float) the angle for the rotation
            part: (None or int>=0) if an integer is provided, only treat
                that part. Please see ManipulateAggregates.aggregate.agg.tag_parts for a
                definition of "part".
        """
        if part is None:
            self.obmol.Rotate(int(axis_index), float(angle))
        else:
            self.obmol.RotatePart(int(part), int(axis_index), float(angle))

    def vdw_check(self, factor=1.0):
        """Check for van-der-Waals clashes.
        
        Check whether any two atoms are closer together than the sum of their
        van-der-Vaals radii.  Perform this check only for atoms that are not
        connected by an arbitrary number of bonds. Hence, this only makes sense
        for aggregates.

        Args:
            factor: (float) before checking for clashes, each
                van-der-Waals radius is multiplied by this factor

        Returns:
            whether or not any two atoms clash
        """
        return self.obmol.IsGoodVDW(float(factor), 0.0)

    def translate(self, vector, part=None):
        """Translate the molecule in a given direction.
    
        Args:
            vector: (list of 3 floats) vector that is added to every atom's coordinate
            part: (None or int>0=) if an integer is provided, only treat that
                part. Please see ManipulateAggregates.aggregate.agg.tag_parts
                for a definition of "part".
        """
        vec = _double_array(vector)
        if part is None:
            self.obmol.Translate(vec)
        else:
            self.obmol.TranslatePart(int(part), vec)
        del vec

    def move_closer(
        self, part1, part2, stepsize=0.1, vdw_factor=0.9, vdw_added=0.0, vec=None
    ):
        """Move two parts of an aggregate closer together.
        
        The parts are moved closer together until a clash of vdW surfaces
        occurs.  Indices start at 1. A "part" is one covalently bound unit as
        determined by the force field. If tagging was enabled using
        ManipulateAggregates.aggregate.agg.tag_parts, a "part" is one
        tagged unit.

        Args:
            part1: (int) index indicating the first part in the aggregate that
                shall be moved closer to another one
            part2: (int) second index
            stepsize: (float) stepsize for movement (good value: 0.2)
            vdw_factor: (float) factor by which all vdW-radii will be
                multiplied before detecting clashes (default: 0.9)
            vdw_added: (float) value that is added to all vdW-radii before
                detecting clashes (default: 0.0)
            vec: (None or list of 3 floats) if not None, the parts will be
                moved closer together in the direction of this vector.
                Otherwise, the parts will be moved closer along the vector
                connecting their non-mass-weighted centers.
        """
        if vec is None:
            self.obmol.MovePartsCloser(
                int(part1),
                int(part2),
                float(stepsize),
                float(vdw_factor),
                float(vdw_added),
            )
        else:
            vtemp = _double_array(vec)
            self.obmol.MovePartsCloser(
                vtemp,
                int(part1),
                int(part2),
                float(stepsize),
                float(vdw_factor),
                float(vdw_added),
            )
            del vtemp

    def tag_parts(self, parts, verbose=True):
        """Manage the use of tagging (and thus change the definition of "part").

        Enable, disable or manage tagging. Part indices start at 0. A "part" is
        one covalently bound unit as determined by the force field. If tagging
        was enabled, a "part" is one tagged unit.

        If parts is a list of integers, the molecules that are indexed with
        the numbers in this list are together added to one tag and tagging is
        enabled for the entire aggregate. If parts is <0, tagging is disabled.
        If parts is >0, tagging is enabled (and nothing else happens). If parts
        is 0, only tagging information is printed.

        Args:
            parts: (int or list of ints) see detailed description
            verbose: (bool) whether or not information about what was changed
                shall be printed

        Returns:
            the index that, from now on, can be used to treat the here tagged
            unit, or None if no tagging was changed.
        """
        try:
            piter = iter(parts)
            iterable = True
        except TypeError:
            iterable = False
        if iterable:
            self.obmol.EnableTags()
            newtag = self.obmol.CreateTag(len(parts))
            for p in piter:
                self.obmol.AddToTag(int(p), newtag)
            if verbose:
                print("Updated tagging information:")
                self.obmol.PrintTags()
        elif parts < 0:
            self.obmol.DisableTags()
            if verbose:
                print("Disabled tagging.")
        elif parts > 0:
            self.obmol.EnableTags()
            if verbose:
                print("Enabled tagging.")
        elif parts == 0:
            self.obmol.PrintTags()

    def append(self, agg, vector=[0, 0, 0], axis=[1, 0, 0], angle=0):
        """Append a molecule to the current one.

        "Appending" means that all the atoms and bonds contained within agg
        are added (deep-dopied, i.e., you can detele agg afterwards) to the
        current molecule. No new bonds are formed. Before appending, translate
        and rotate agg.

        Args:
            agg: aggregate to be appended (of the same type)
            vector: (list of 3 floats) translate agg by this vector prior to appending
            axis: (list of 3 floats) rotate agg around this axis prior to appending it
            angle: (float) the angle for the rotation

        """
        vec = _double_array(vector)
        ax = _double_array(axis)
        self.obmol.AppendMolecule(agg.obmol, vec, ax, float(angle))
        del vec, ax

    def glue(self, agg, i1, i2, m1, m2):
        """Glue a molecule to the current one.
        
        Two bonds, one in each involved molecule, will be cut in half and then
        glued together in such a way, that all atoms connected to i2 and
        m2 (including those two) but not connected to i1 and m1 will
        be cleaved and the molecules will be glued together so that i1 and
        m1 form a proper bond.
    
        Args:
            agg: aggregate to be glued to this one (of same type)
            i1: (int) index specifying the atom of the primary molecule that
                will be retained
            i2: (int) index specifying the atom of the primary molecule that
                will not be retained
            m1: (int) index specifying the atom of the secondary molecule that
                will be retained
            m2: (int) index specifying the atom of the secondary molecule that
                will not be retained
        """
        self.obmol.BondMolecule(agg.obmol, int(i1), int(i2), int(m1), int(m2))

    def cleave(self, i1, i2):
        """Cleave a part of a molecule.

        Cleave all atoms and bonds that are connected to the atom indexed by
        i2 but not to i1. Leave the rest as it is.
    
        Args:
            i1: (int) index specifying the atom of the molecule that
                will be retained
            i2: (int) index specifying the atom of the molecule that
                will not be retained
        """
        self.obmol.CleaveOff(int(i1), int(i2))

    def write(self, filename, fileformat=None, overwrite=True):
        """Write the data of the molecule to disk. 
        
        Args:
            filename: (string) path to the file (INCLUDING the extension)
            fileformat: (string) output file format (anything that maagbel can write)
            overwrite: shall the output file be overwritten or not

        Raises:
            ManipulateAggregates.aggregate.MissingModuleError.
        """
        if fileformat is None:
            fileformat = guess_format(filename, True)
        if fileformat is None:
            fileformat = self.info.get("outformat", None)
        if fileformat is None:
            raise FiletypeException(
                "Filetype of %s could not be automatically determined and this aggregate was not created from a file."
                % (filename)
            )
        pybel.Molecule(self.obmol).write(fileformat, filename, overwrite=overwrite)

    def align(self, point, main3, main2, part=None):
        """Align an aggregate in space.

        Align the last two main axes of an aggregate to the two given axes and
        move the center to the given coordinate.
    
        Args:
            point: (list of 3 floats) the new center of the molecule (not mass weighed)
            main3: (list of 3 floats) the new 3rd main axis (usually the longest extent)
            main2: (list of 3 floats) the new 2nd main axis
            part: (None or int>=0) if an integer is provided, only treat that
                part. Please see ManipulateAggregates.aggregate.agg.tag_parts
                for a definition of "part".
        """
        for vec in [[point, "point"], [main3, "third axis"], [main2, "second axis"]]:
            if not len(vec[0]) == 3:
                raise IndexError(
                    "Variable "
                    + vec[1]
                    + " not of the correct length, needs 3 elements not "
                    + str(len(vec[0]))
                )
        if sum(abs(v) for v in main3) > 0 and sum(abs(v) for v in main2) > 0:
            poi = _double_array(point)
            ma3 = _double_array(main3)
            ma2 = _double_array(main2)
            if part is None:
                self.obmol.Align(poi, ma3, ma2)
            else:
                self.obmol.AlignPart(int(part), poi, ma3, ma2)
            del poi, ma3, ma2
        else:
            raise ValueError("Main axis vectors have to have a non-vanishing norm.")

    def mirror(self, normal, point, center_it=False, part=None):
        """Mirror or point-invert an aggregate.

        Mirror the molecule either by point inversion or by mirroring at a
        plane.
    
        Args:
            normal: (list of 3 floats) the normal vector of the mirror plane.
                If its norm is 0, point inversion will be performed.
            point: (list of 3 floats) either the inversion point or a point in
                the plane (Hessian normal form)
            center_it: (bool) whether or not the entire aggregate's center
                shall be moved to the origin prior to the operation. The center
                will automatically be moved back to its original position
                afterward the procedure.
            part: (None or int>=0) if an integer is provided, only treat that
                part. Please see ManipulateAggregates.aggregate.agg.tag_parts
                for a definition of "part".
        """
        nor = _double_array(normal)
        poi = _double_array(point)
        if part is None:
            self.obmol.Mirror(nor, poi, center_it)
        else:
            self.obmol.MirrorPart(int(part), nor, poi, center_it)
        del nor, poi

    def part(self, normal_vector, coordinate):
        """Get an OBMol containing all those atoms that are one one side of a plane.

        The plane is given in the Hessian normal form.

        Args:
            normal_vector: (list of 3 floats) the normal vector of the plane
            coordinate: (list of 3 floats) a point in the plane
        """
        tempmol = maagbel.OBMol()
        nor = _double_array(normal_vector)
        coo = _double_array(coordinate)
        self.obmol.PartMolecule(tempmol, nor, coo)
        del nor, coo
        return tempmol

    def part_aggregate(self, normal_vector, coordinate, side="left"):
        """
        Get an object of the same type containing all those atoms that are one one side of a plane.

        The plane is given in the Hessian normal form. Please note that the
        contructor will be called.

        Args:
            normal_vector: (list of 3 floats) the normal vector of the plane
            coordinate: (list of 3 floats) a point in the plane
            side: (string) If 'left', select those atoms on the side where the
                normal vector points. If anything else, select all those on the
                opposite side.
        """
        if side == "right":
            normal_vector = [-i for i in normal_vector]
        elif side == "left":
            pass
        else:
            raise ValueError("Side must be either left or right")
        return agg(self.part(normal_vector, coordinate), ff=None)

    def get_partial_charges(self):
        """Get a list of all partial charges of the atoms according to the specified method.
        
        The returned list has the same order as the atoms in the molecule.

        Returns:
            a list of floats containing the partial charges.

        Raises:
            ValueError, ManipulateAggregates.aggregate.OpenBabelError.
        """
        method = self.cp["method"].lower()
        if method == "corecharge":
            partialcharges = self.get_charges()
        else:
            tmp_charges = maagbel.OBChargeModel.FindType(method)
            if tmp_charges is not None:
                if tmp_charges.ComputeCharges(self.obmol):
                    partialcharges = list(tmp_charges.GetPartialCharges())
                else:
                    raise OpenBabelError("Error while partitioning partial charges.")
            else:
                raise ValueError(
                    "Method '"
                    + method
                    + "' is not a known method for partitioning partial charges. See 'manipagg --list charges' for partitioning methods or use 'corecharge' to use the core charge."
                )
            del tmp_charges
        # qeq does deliver charges of the opposite sign as the rest
        if method == "qeq":
            partialcharges = [-q for q in partialcharges]
        return partialcharges

    def set_ecp(self, idx, ecp):
        """Set an atom's (or of multiple ones) core charge property.

        If idx and ecp are iterables, idx-ecp pairs will be assigned.
        If the iterables are of unequal lenghts (e.g., n and n+m), only the
        first n key-value pairs will be treated.

        Args:
            idx: (int or iterable of ints) the id(s) (starting at 1) of the
                atom whose ecp shall be set
            ecp: (int or iterable of ints) the cored charge(s) of the atom
                whose ecp shall be set
        """
        # a=maagbel.OBAtom()
        try:
            ieiter = iter(zip(idx, ecp))
            iterable = True
        except TypeError:
            iterable = False
        if iterable:
            for i, e in ieiter:
                if e != 0:
                    a = self.obmol.GetAtom(int(i))
                    pd = maagbel.OBPairData()
                    pd.SetAttribute("ecp")
                    pd.SetValue(str(e))
                    a.CloneData(pd)
        else:
            if ecp != 0:
                a = self.obmol.GetAtom(int(idx))
                pd = maagbel.OBPairData()
                pd.SetAttribute("ecp")
                pd.SetValue(str(ecp))
                a.CloneData(pd)

    def get_charges(self, ecp=True):
        """Return a list of all charges of the atoms according to ther element numbers.

        This function respects the atom property "ecp" that might be stored in
        each of the OBAtom objects, if ecp is True.

        Args:
            ecp: (bool) whether or not to respect the core charge property
        
        Returns:
            a list of floats containing the elemental charges (possibly minus core charges).
        """
        # a=maagbel.OBAtom()
        charges = [0.0] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            charges[idx - 1] = a.GetAtomicNum()
            if ecp:
                aecp = a.GetData("ecp")
                if aecp is not None:
                    charges[idx - 1] -= int(aecp.GetValue())
        return charges

    def get_dipole_moment(self):
        """Get an aggregate's electric dipole moment from point charges.

        The currently set method for obtaining partial charges will be used.

        Returns:
            a list of 3 floats, the electric dipole vector
        """
        charges = self.get_partial_charges()
        coordinates = self.get_coordinates()
        px = 0
        py = 0
        pz = 0
        for (x, y, z), c in zip(coordinates, charges):
            px += x * c
            py += y * c
            pz += z * c
        del charges, coordinates
        return [px, py, pz]

    def get_coordinates(self):
        """Get the coordinates of all atoms.

        Returns:
            a list of lists of 3 floats, the Cartesian coordinates of the atoms
        """
        # a=maagbel.OBAtom()
        coordinates = [None] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            coordinates[idx - 1] = [a.GetX(), a.GetY(), a.GetZ()]
        return coordinates

    def get_center(self, mask=None):
        """Get the non-mass-weighted center of the molecule

        Returns:
            a list of 3 floats, the Cartesian coordinates of the center
        """
        _assert_supported("numpy")
        coords = self.get_coordinates()
        if mask is not None:
            if True in (i <= 0 for i in mask):
                raise ValueError("Counting atoms starts at 1, even for masks.")
            coords = [coords[i - 1] for i in mask]
        return numpy.mean(numpy.array(coords), axis=0)

    def get_obatom_vec(self, mask=None):
        """Get a vector of OBAtom pointers

        Returns:
            a swig proxy to the vector of atoms
        """
        atomlist = []
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            atomlist.append(a)
        if mask is not None:
            if True in (i <= 0 for i in mask):
                raise ValueError("Counting atoms starts at 1, even for masks.")
            atomlist = [atomlist[i - 1] for i in mask]
        obatoms = maagbel.vectorOBAtom(atomlist)
        return obatoms

    def get_main_axes(self, mask=None):
        """Get the last 2 main axes of the aggregate.

        Returns:
            a tuple of 2 lists of 3 floats, the third and second main axes
        """
        _assert_supported("numpy")
        # center = numpy.array(self.get_center(mask))
        # coords = numpy.array(self.get_coordinates())-center
        if mask is not None:
            if True in (i <= 0 for i in mask):
                raise ValueError("Counting atoms starts at 1, even for masks.")
            bvmask = self.get_bit_mask(mask)
        ##mat is the tensor of inertia
        # mat    = numpy.sum(numpy.array([ [[y*y+z*z,-x*y,-x*z],[-x*y,x*x+z*z,-y*z],[-x*z,-y*z,x*x+y*y]] for x,y,z in coords]),axis=0)
        # eigvals,eigvecs = numpy.linalg.eig(mat)
        ##the eigenvectors are stored in the coloumns of eigvecs
        ##so it is transposed to have easy access to them
        # eigvecs = -eigvecs.T
        # t_main3,t_main2,t_main1 = sorted(zip(eigvals,eigvecs),key=lambda e: e[0])

        p = maagbel.vector3(0, 0, 0)
        m3 = maagbel.vector3(0, 0, 0)
        m2 = maagbel.vector3(0, 0, 0)
        if mask is not None:
            self.obmol.GetMainAxes(p, m3, m2, bvmask)
        else:
            self.obmol.GetMainAxes(p, m3, m2)
        main3 = (m3.GetX(), m3.GetY(), m3.GetZ())
        main2 = (m2.GetX(), m2.GetY(), m2.GetZ())

        return main3, main2

    def get_align_matrix(self, main3, main2):
        """Get a matrix to align the aggregate.

        Return the composite rotation matrix that would align the third and
        second main axes to the given axes.

        Args:
            main3: (list of 3 floats) the new 3rd main axis (usually the longest extent)
            main2: (list of 3 floats) the new 2nd main axis

        Returns:
            a numpy array of shape (3,3) and dtype float, the rotation matrix
        """
        _assert_supported("numpy")
        # c_ stands for current
        c_main3, c_main2 = self.get_main_axes()

        tempvec = numpy.cross(main3, c_main3)
        angle = _VectorAngle(main3, c_main3)
        mat1 = _RotMatrixAboutAxisByAngle(tempvec, angle)

        c_main2 = numpy.dot(mat1, c_main2)
        tempvec = numpy.cross(main2, c_main2)
        angle = _VectorAngle(main2, c_main2)
        mat2 = _RotMatrixAboutAxisByAngle(tempvec, angle)

        return numpy.dot(mat2, mat1)

    def get_povlight_matrix(self):
        """Get a matrix using which to rotate all normal vectors prior to povray visualization.

        Return the composite rotation matrix that allow a rotation around an
        axis by an angle. This depends on the currect configuration option for
        the key "visrotmat".

        Returns:
            a numpy array of shape (3,3) and dtype float, the rotation matrix
        """
        _assert_supported("numpy")
        axis, angle = self.vs.get("visrotmat", ([1, 0, 0], 0))
        mat = _RotMatrixAboutAxisByAngle(axis, numpy.pi * angle / 180)
        return mat

    def get_vdw_radii(self):
        """Get all van-der-Waals radii for the atoms in this aggregate.

        Returns:
            a list of floats, van-der-Waals radii of the atoms according to the
            element table of OpenBabel
        """
        # a=maagbel.OBAtom()
        vdw_radii = [0.0] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            vdw_radii[idx - 1] = maagbel.etab.GetVdwRad(a.GetAtomicNum())
        return vdw_radii

    def get_names(self):
        """Get a list of all element symbols of the atoms.

        Returns:
            a list of strings, the element symbols of the atoms.
        """
        # a=maagbel.OBAtom()
        names = [None] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            names[idx - 1] = maagbel.etab.GetSymbol(a.GetAtomicNum())
        return names

    def get_colours(self):
        """Get a list of colors associated with the element symbols of the atoms.

        Returns:
            a list of tuples of 3 floats, RGB values
        """
        # a=maagbel.OBAtom()
        colours = [None] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            colours[idx - 1] = maagbel.etab.GetRGB(a.GetAtomicNum())
        return colours

    def get_masses(self):
        """Get a list of masses associated with the element symbols of the atoms.

        Returns:
            a list of floats, the masses in atomic units
        """
        # a=maagbel.OBAtom()
        masses = [0.0] * self.obmol.NumAtoms()
        for idx in range(1, self.obmol.NumAtoms() + 1):
            a = self.obmol.GetAtom(int(idx))
            masses[idx - 1] = maagbel.etab.GetMass(a.GetAtomicNum())
        return masses

    def get_vdw_surface(self, nr_refinements=1, shrink_factor=0.95, vdwscale=1.0):
        """Compute the aggregate's discretized van-der-Waals surface.

        Do not call directly. Rather use ManipulateAggregates.aggregate.agg.get_surface
        and set properties via ManipulateAggregates.aggregate.agg.set_vs beforehand.

        Uses CGAL's skin surface mesher. Please see
        http://doc.cgal.org/latest/Skin_surface_3/index.html for more
        information on this.

        Args:
            nr_refinements: (int) number of refinement steps for the skin
                surface generation. The higher the number the more vertices it
                will have.
            shrink_factor: the shrink factor for the generation of the skin
                surface. Must be between 0 and 1 (not including those). The
                bigger the value the tighter the surface will be.
            vdwscale: (float) multiply each vdW radius by this value before
                building the surface

        Returns:
            a tuple of corners,face_indices,normals. corners (a list of lists
            of 3 floats) contains the Cartesian coordinates of all vertices of
            the surface. face_indices (list of lists of 3 ints) each triple of
            integers defines one face of the surface. The indices correspond to
            corners. normals (list of lists of 3 floats) each triple defines
            the normal vector associated with the corresponding face.


        Example in 2D for nr_points=12

        . : point on the sphere's surface
        
        X : center of the sphere

        Not overlapping => 12 points per sphere

        >>>    ...   ...
        >>>   .   . .   .
        >>>   . X . . X .
        >>>   .   . .   .
        >>>    ...   ...

        Overlapping => points that would be within the other sphere are removed => 9
        points per "sphere"

        >>>    ......
        >>>   .      .
        >>>   . X  X .
        >>>   .      .
        >>>    ......

        """
        _assert_supported("numpy")
        _assert_supported("FireDeamon")

        vdw_radii = [r * vdwscale for r in self.get_vdw_radii()]
        coordinates = self.get_coordinates()

        lengths, face_indices, corners, normals = fd.SkinSurfacePy(
            shrink_factor, coordinates, vdw_radii, refinesteps=nr_refinements
        )

        return corners, face_indices, normals

    def get_iso_surface(
        self,
        isovalue,
        isofile,
        isofiletype="dx",
        mesh_criteria=MESH_CRITERIA,
        relative_precision=1.0e-06,
        atoms="auto",
    ):
        """Compute a discretized iso surface of the aggregate.

        Do not call directly. Rather use ManipulateAggregates.aggregate.agg.get_surface
        and set properties via ManipulateAggregates.aggregate.agg.set_vs beforehand.
        Please see http://doc.cgal.org/latest/Surface_mesher/index.html for
        more information.

        Args:
            isovalue: (float) the isovalue for the surface
            isofile: (string) path to the file from which to take the volumetric data

        Args:
            isofiletype: (string) filetype of the volumetric data file. Only "dx" is
                supported as of now. "cube" might be added later.
            mesh_criteria: (list of 3 floats) CGAL's internal meshing criteria
            relative_precision: (float) CGAL's internal precision for meshing
            atoms: (int, list of ints or "all" or "noH" or "auto") CGAL's mesh
                generation requires a point inside the isosurface. For a
                chemical compound, all atoms should lie within the isosurface.
                For aggregates, this is not the case. Here, specofy at least
                one atom of every covalently bound unit. The special values
                "all" and "noH" select all atoms or only non-hydrogen atoms.
                The special value "auto" automatically uses the first atom of
                each covalently bound unit.

        Returns:
            a tuple of corners,face_indices,normals. corners (a list of lists
            of 3 floats) contains the Cartesian coordinates of all vertices of
            the surface. face_indices (list of lists of 3 ints) each triple of
            integers defines one face of the surface. The indices correspond to
            corners. normals (list of lists of 3 floats) each triple defines
            the normal vector associated with the corresponding face.
        """
        _assert_supported("numpy")
        _assert_supported("FireDeamon")
        _assert_supported("ManipulateAggregates.collection.read")

        if atoms == "auto":
            # GetConnections gives back a string that contains
            # all the atom indices in one covalently bound unit
            # per line. The first entry, however, is the molecule
            # index.
            s = self.obmol.GetConnections()
            # split by line -> convert each entry to int (apart form first)
            #               -> sort by number -> take first entry of sorted list
            atoms = [
                m[0]
                for m in (
                    sorted(map(int, e.split()[1:])) for e in s.split("\n") if len(e) > 0
                )
            ]
            if len(atoms) == 1:
                atoms = atoms[0]

        if isinstance(atoms, int):
            print(
                "WARNING: Using only one atom to generate iso surface.", file=sys.stderr
            )
            coordinates = [self.get_coordinates()[atoms]]
        elif atoms == "noH":
            coordinates = [
                c
                for c, n in zip(self.get_coordinates(), self.get_names())
                if not n == "H"
            ]
        elif atoms == "all":
            coordinates = self.get_coordinates()
        else:
            coordinates = self.get_coordinates()
            coordinates = [coordinates[i] for i in atoms]

        if isofiletype.lower() == "dx":
            header = {}
            data = fdread.read_dx(
                isofile,
                unit_conversion=1.0,
                invert_charge_data=False,
                density=True,
                header_dict=header,
                grid=False,
                data=True,
                silent=False,
                gzipped=False,
                comments=False,
            )
            origin = header["org_xyz"]
            counts = header["counts_xyz"]
            delta = [header["delta_x"], header["delta_y"], header["delta_z"]]
            data = data["data"]

            lengths, face_indices, corners, normals = fd.IsosurfacePy(
                data,
                origin,
                counts,
                delta,
                isovalue,
                coordinates,
                relative_precision,
                mesh_criteria,
            )
        else:
            raise ValueError(
                "Type of isofile '%s' unknown. Supported are: 'dx'." % (isofiletype)
            )

        return corners, face_indices, normals

    def get_surface(self):
        """Compute a discretized surface of the aggregate according to the current config.

        Set properties via ManipulateAggregates.aggregate.agg.set_vs beforehand. See
        ManipulateAggregates.aggregate.agg.default_vs for config options.

        Returns:
            a tuple of corners,face_indices,normals. corners (a list of lists
            of 3 floats) contains the Cartesian coordinates of all vertices of
            the surface. face_indices (list of lists of 3 ints) each triple of
            integers defines one face of the surface. The indices correspond to
            corners. normals (list of lists of 3 floats) each triple defines
            the normal vector associated with the corresponding face.
        """
        if self.vs["type"].lower() == "vdw":
            return self.get_vdw_surface(
                nr_refinements=self.vs["refine"],
                shrink_factor=self.vs["shrink_factor"],
                vdwscale=self.vs["vdw_scale"],
            )
        elif self.vs["type"].lower() == "iso":
            return self.get_iso_surface(
                self.vs["isovalue"],
                self.vs["isofile"],
                isofiletype=self.vs["isofiletype"],
                mesh_criteria=self.vs["mesh_criteria"],
                relative_precision=self.vs["rel_precision"],
                atoms=self.vs["iso_atoms"],
            )
        else:
            raise ValueError("Only iso and van-der-Waals surfaces are supported.")

    def get_density(self, points):
        """Compute this aggregate's electron density at the given coordinates.

        This function checks whether the last call used the same configuration
        as the previous one. If that is so, files are not read in again.

        Args:
            points: (list of lists of 3 floats) the Cartesian coordinates at
                which to compute the density 

        Returns:
            a numpy array of dtype float containing the density values at the
            specified points
        """

        prog_report = os.environ.get("PROGRESS", True)
        if prog_report == "0":
            prog_report = False
        else:
            prog_report = True

        _assert_supported("numpy")
        _assert_supported("FireDeamon")

        corners = numpy.array(points)

        ################################################################################
        if self.cp["type"] == "empirical":
            raise ValueError(
                "Cannot compute electron density due to empirical charges."
            )
        ################################################################################
        elif self.cp["type"] == "orbitals":
            orbcfg = {
                "orbfiletype": self.cp["orbfiletype"].lower(),
                "orbfile": self.cp["orbfile"],
            }
            refresh = not (self.__internal__["orbcfg"] == orbcfg)
            if refresh:
                self.__internal__["orbcfg"] = orbcfg
            ###########
            if self.cp["orbfiletype"].lower() == "molden":
                BOHRTOANG = orbitalcharacter.BOHRTOANG
                if refresh:
                    self.__internal__[
                        "orb"
                    ] = orbitalcharacter.read_molden_orbitals_corrected(
                        self.cp["orbfile"]
                    )
                (
                    basis,
                    Smat,
                    (MOsalpha, MOsbeta),
                    (OCCsalpha, OCCsbeta),
                ) = self.__internal__["orb"]
                data = fd.InitializeGridCalculationOrbitalsPy(
                    corners, basis, scale=BOHRTOANG
                )
            else:
                raise ValueError("Unkown orbital file type.")
            ###########
            if MOsalpha == MOsbeta:
                density = numpy.array(
                    fd.ElectronDensityPy(
                        MOsalpha,
                        data,
                        occupations=[2 * o for o in OCCsalpha],
                        prog_report=prog_report,
                        cutoff=self.cp["cutoff"],
                    )
                )
            else:
                density = numpy.array(
                    fd.ElectronDensityPy(
                        MOsalpha + MOsbeta,
                        data,
                        occupations=OCCsalpha + OCCsbeta,
                        prog_report=prog_report,
                        cutoff=self.cp["cutoff"],
                    )
                )
        ################################################################################
        elif self.cp["type"] == "interpolation":
            raise ValueError("Cannot interpolate electron density.")
        ################################################################################
        elif self.cp["type"] == "charges":
            raise ValueError("Cannot compute electron density from charges.")
        ################################################################################
        else:
            raise ValueError(
                "Unknown value for key 'type' %s for obtaining the electron density. I know: 'orbitals'."
                % (self.cp["type"])
            )

        return density

    def get_potential(self, points):
        """Compute this aggregate's electrostatic potential at the given coordinates.

        This function checks whether the last call used the same configuration
        as the previous one. If that is so, files are not read in again.

        Args:
            points: (list of lists of 3 floats) the Cartesian coordinates at
                which to compute the potential

        Returns:
            a numpy array of dtype float containing the potential values at the
            specified points
        """

        prog_report = os.environ.get("PROGRESS", True)
        if prog_report == "0":
            prog_report = False
        else:
            prog_report = True

        _assert_supported("numpy")
        _assert_supported("FireDeamon")

        corners = numpy.array(points)

        ################################################################################
        if self.cp["type"] == "empirical":
            charges = self.get_partial_charges()
            coordinates = self.get_coordinates()
            if not self.cp["partial"]:
                charges = [c + cc for c, cc in zip(charges, self.get_charges())]
            potential = numpy.array(
                fd.ElectrostaticPotentialPy(
                    points, charges, coordinates, prog_report=prog_report
                )
            )
        ################################################################################
        elif self.cp["type"] == "orbitals":
            orbcfg = {
                "orbfiletype": self.cp["orbfiletype"].lower(),
                "orbfile": self.cp["orbfile"],
            }
            refresh = not (self.__internal__["orbcfg"] == orbcfg)
            if refresh:
                self.__internal__["orbcfg"] = orbcfg
            ###########
            if self.cp["orbfiletype"].lower() == "molden":
                BOHRTOANG = orbitalcharacter.BOHRTOANG
                if refresh:
                    self.__internal__[
                        "orb"
                    ] = orbitalcharacter.read_molden_orbitals_corrected(
                        self.cp["orbfile"]
                    )
                (
                    basis,
                    Smat,
                    (MOsalpha, MOsbeta),
                    (OCCsalpha, OCCsbeta),
                ) = self.__internal__["orb"]
                data = fd.InitializeGridCalculationOrbitalsPy(
                    corners, basis, scale=BOHRTOANG
                )
            else:
                raise ValueError("Unkown orbital file type.")
            ###########
            if MOsalpha == MOsbeta:
                potential = -numpy.array(
                    fd.ElectrostaticPotentialOrbitalsPy(
                        MOsalpha,
                        Smat,
                        [2 * o for o in OCCsalpha],
                        data,
                        prog_report=prog_report,
                    )
                )
            else:
                potential = -numpy.array(
                    fd.ElectrostaticPotentialOrbitalsPy(
                        MOsalpha + MOsbeta,
                        Smat,
                        OCCsalpha + OCCsbeta,
                        data,
                        prog_report=prog_report,
                    )
                )

            charges = self.get_charges()
            coordinates = self.get_coordinates()
            pospotential = numpy.array(
                fd.ElectrostaticPotentialPy(
                    corners / BOHRTOANG,
                    charges,
                    [[xyz / BOHRTOANG for xyz in a] for a in coordinates],
                    prog_report=prog_report,
                )
            )
            potential += pospotential
        ################################################################################
        elif self.cp["type"] == "interpolation":
            potcfg = {
                "potfiletype": self.cp["potfiletype"].lower(),
                "potfile": self.cp["potfile"],
                "total_charge": self.cp["total_charge"],
                "method": self.cp["interpolation"].lower(),
                "function": self.cp["int_root"],
                "exponent": self.cp["int_exponent"],
                "cutoff": self.cp["cutoff"],
            }
            refresh = not (self.__internal__["potcfg"] == potcfg)
            if refresh:
                self.__internal__["potcfg"] = potcfg
                ###########
                if self.cp["potfiletype"].lower() == "xyz":
                    coordinates, potential = fdread.read_charges_simple(
                        self.cp["potfile"]
                    )
                ###########
                elif self.cp["potfiletype"].lower() == "dx":
                    coordinates, potential = fdread.read_charges_dx(
                        self.cp["potfile"],
                        add_nuclear_charges=False,
                        unit_conversion=1.0,
                        total_charge=self.cp["total_charge"],
                        invert_charge_data=False,
                        rescale_charges=False,
                        density=True,
                    )
                ###########
                elif self.cp["potfiletype"].lower() == "cube":
                    coordinates, charges = read_charges_cube(
                        self.cp["potfile"],
                        add_nuclear_charges=False,
                        rescale_charges=False,
                        density=True,
                    )
                ###########
                else:
                    raise ValueError("Unkown potential file type.")
                ##########
                self.__internal__["pot"] = (coordinates, potential)
            else:
                coordinates, charges = self.__internal__["pot"]
            ###########
            config = {
                "method": self.cp["interpolation"].lower(),
                "function": self.cp["int_root"],
                "exponent": self.cp["int_exponent"],
                "cutoff": self.cp["cutoff"],
            }
            potential = numpy.array(
                fd.InterpolationPy(
                    coordinates,
                    potential,
                    points,
                    prog_report=prog_report,
                    config=config,
                )
            )
        ################################################################################
        elif self.cp["type"] == "charges":
            chacfg = {
                "chargefiletype": self.cp["chargefiletype"].lower(),
                "chargefile": self.cp["chargefile"],
                "partial": self.cp["partial"],
                "total_charge": self.cp["total_charge"],
            }
            refresh = not (self.__internal__["chacfg"] == chacfg)
            if refresh:
                self.__internal__["chacfg"] = chacfg
                ###########
                if self.cp["chargefiletype"].lower() == "xyz":
                    if self.cp["partial"]:
                        coordinates, charges = fdread.read_charges_simple(
                            self.cp["chargefile"], compare_elements=True, molecule=self
                        )
                    else:
                        coordinates, charges = fdread.read_charges_simple(
                            self.cp["chargefile"]
                        )
                ###########
                elif self.cp["chargefiletype"].lower() == "dx":
                    coordinates, charges = fdread.read_charges_dx(
                        self.cp["chargefile"],
                        add_nuclear_charges=True,
                        molecule=self,
                        unit_conversion=1.0,
                        total_charge=self.cp["total_charge"],
                        invert_charge_data=True,
                    )
                ###########
                elif self.cp["chargefiletype"].lower() == "cube":
                    coordinates, charges = read_charges_cube(
                        self.cp["chargefile"],
                        add_nuclear_charges=True,
                        total_charge=self.cp["total_charge"],
                    )
                ###########
                else:
                    raise ValueError("Unkown charge file type.")
                ###########
                self.__internal__["cha"] = (coordinates, charges)
            else:
                coordinates, charges = self.__internal__["cha"]
            ###########
            potential = numpy.array(
                fd.ElectrostaticPotentialPy(
                    points, charges, coordinates, prog_report=prog_report
                )
            )
        ################################################################################
        else:
            raise ValueError(
                "Unknown value for key 'type' %s for obtaining the electrostatic potential. I know: 'empirical', 'orbitals', 'interpolation' and 'charges'."
                % (self.cp["type"])
            )

        if self.cp["invert_potential"]:
            potential *= -1.0

        return potential

    def get_bond_map(self, unique=True, no_hydrogen=False):
        """Produce a list of all bonds in a molecule as known by the current force field.

        Returns:
            a list of tuples fo 2 ints, the atom indices associated with the bonds

        Args:
            unique: (bool) if True, give back an irreducible list of bonds in the form
                of tuples of indices.  If False, give back a complete list of bonds, i.e.
                every atom in a bond is once the first and once the second element in one of
                the tuples

            no_hydrogen: (bool) whether or not to exclude hydrogens from the list
        """
        bondmap = []
        for bond_id in range(0, self.obmol.NumBonds()):
            b = self.obmol.GetBond(int(bond_id))
            if (
                no_hydrogen
                and b.GetBeginAtom().IsHydrogen()
                or b.GetEndAtom().IsHydrogen()
            ):
                continue
            bondmap.append((b.GetBeginAtomIdx() - 1, b.GetEndAtomIdx() - 1))
        if unique:
            bondmap = [tuple(sorted(b, key=lambda x: x)) for b in bondmap]
            bondmap = sorted(
                list(set(bondmap)), key=lambda x: x[0] * (len(bondmap) + 1) + x[1]
            )
        else:
            bondmap = sorted(bondmap, key=lambda x: x[0] * (len(bondmap) + 1) + x[1])
        return bondmap

    def visualize(self):
        """Visualize the aggregate according to the currect config."""
        visualize.visualize(self)

    def rmsd(self, agg, print_result=False):
        """Determine the Root-Mean-Square-Deviation with respect to another aggregate.

        Args:
            agg: (of same type) the aggregate to compare agains
            print_result: (bool) whether or not to also print the result

        Raises:
            ValueError.
        """
        _assert_supported("numpy")
        selftemp = self.duplicate()
        othertemp = agg.duplicate()
        selftemp.align([0, 0, 0], [1, 0, 0], [0, 1, 0])
        othertemp.align([0, 0, 0], [1, 0, 0], [0, 1, 0])
        selfcoords = numpy.array(selftemp.get_coordinates())
        othercoords = numpy.array(othertemp.get_coordinates())
        if not (selfcoords.shape == othercoords.shape):
            raise ValueError(
                "The aggregate to compare against does not have the same number of atoms."
            )
        diff = selfcoords - othercoords
        result_rmsd = numpy.sqrt(numpy.sum(diff * diff) / (3.0 * len(selfcoords)))
        result_maxdeviation_single = numpy.max(numpy.abs(diff))
        result_maxdeviation_whole = numpy.max(numpy.linalg.norm(diff, axis=1))
        del selfcoords, othercoords
        del selftemp, othertemp
        result = (result_rmsd, result_maxdeviation_single, result_maxdeviation_whole)
        if print_result:
            print(
                "RMSD: %.2e | Max deviation in one coordinate: %.2e | Max deviation for a single atom: %.2e"
                % result
            )
        return result
