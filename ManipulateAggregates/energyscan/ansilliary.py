"""A set of ansilliary functions for the energyscan module.
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
import sys
import os
import re
import errno

import logging

logger = logging.getLogger(__name__)

try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from maagbel import doubleArray
except ImportError:
    logger.warning("Could not import doubleArray from maagbel")
try:
    from ..collection.write import print_dx_file as c_print_dx_file
except ImportError:
    logger.warning(
        "Could not import ..collection.write.print_dx_file as c_print_dx_file"
    )
try:
    from ..collection import hashIO
except ImportError:
    logger.warning("Could not import ..collection.read.hashIO")
try:
    FileExistsError("")
except NameError:

    class FileExistsError(OSError):
        pass


# escape sequence to move the cursor up one line
CURSOR_UP_ONE = "\x1b[1A"

# escape sequence to erase the current line
ERASE_LINE = "\x1b[2K"

# conversion factors from meV to the declared force field units
E_UNIT_CONVERSION = {
    "kJ/mol": 0.09648500,
    "kcal/mol": 0.02306035,
}

# associate the name of each pointgroup with its subgroups
# C1 is implicit
SUBGROUPS = {
    "C4v": ("C2v", "C2", "C4", "Cs"),
    "D4h": (
        "C4v",
        "C2v",
        "D2h",
        "C2h",
        "C4h",
        "D2d",
        "C2",
        "C4",
        "S4",
        "Cs",
        "Ci",
        "D4",
        "D2",
    ),
    "C2v": ("C2", "Cs"),
    "D6h": (
        "C2v",
        "D3h",
        "C6v",
        "S6",
        "D2h",
        "C6h",
        "C2h",
        "C3",
        "C2",
        "C6",
        "C3h",
        "Cs",
        "D3d",
        "C3v",
        "Ci",
        "D6",
        "D2",
        "D3",
    ),
    "C6v": ("C2v", "C3", "C2", "C6", "Cs", "C3v"),
    "D6d": ("C2v", "C6v", "D2d", "C3", "C2", "C6", "S4", "Cs", "C3v", "D6", "D2", "D3"),
    "S6": ("C3", "Ci"),
    "D4d": ("C4v", "C2v", "C2", "C4", "Cs", "S8", "D4", "D2"),
    "D2h": ("C2v", "C2h", "C2", "Cs", "Ci", "D2"),
    "S8": ("C2", "C4"),
    "C6h": ("S6", "C2h", "C3", "C2", "C6", "C3h", "Cs", "Ci"),
    "C2h": ("C2", "Cs", "Ci"),
    "C4h": ("C2h", "C2", "C4", "S4", "Cs", "Ci"),
    "D2d": ("C2v", "C2", "S4", "Cs", "D2"),
    "D3d": ("C2v", "S6", "D2h", "C2h", "C3", "C2", "Cs", "C3v", "Ci", "D2", "D3"),
    "C3v": ("C3", "Cs"),
    "C8v": ("C4v", "C2v", "C8", "C2", "C4", "Cs"),
    "C8": ("C2", "C4"),
    "O": ("C3", "C2", "C4", "T", "D4", "D2", "D3"),
    "D8h": (
        "C4v",
        "C8h",
        "C2v",
        "D4d",
        "D2h",
        "D4h",
        "C2h",
        "C4h",
        "D2d",
        "C8v",
        "C8",
        "C2",
        "C4",
        "S4",
        "Cs",
        "Ci",
        "S8",
        "D8",
        "D4",
        "D2",
    ),
    "C8h": ("C2h", "C4h", "C8", "C2", "C4", "S4", "Cs", "Ci", "S8"),
    "C2": (),
    "D8d": ("C4v", "C2v", "C8v", "C8", "C2", "C4", "Cs", "D8", "D4", "D2"),
    "C7": (),
    "C6": ("C3", "C2"),
    "C5": (),
    "C4": ("C2",),
    "D5h": (),
    "C5v": ("C5", "Cs"),
    "C3h": ("C3", "Cs"),
    "Oh": (
        "C4v",
        "C2v",
        "S6",
        "D2h",
        "D4h",
        "C2h",
        "C4h",
        "D2d",
        "C3",
        "C2",
        "C4",
        "Td",
        "S4",
        "Cs",
        "T",
        "D3d",
        "C3v",
        "O",
        "Th",
        "Ci",
        "D4",
        "D2",
        "D3",
    ),
    "I": ("C3", "C2", "C5", "T", "D5", "D2", "D3"),
    "K": (),
    "Kh": ("K", "Cs", "Cinfv", "Ci"),
    "D5d": ("C2v", "D2h", "C2h", "C2", "C5", "C5v", "Cs", "Ci", "D5", "D2"),
    "Td": ("C2v", "D2d", "C3", "C2", "S4", "Cs", "T", "C3v", "D2", "D3"),
    "C7v": ("C7", "Cs"),
    "C3": (),
    "T": ("C3", "C2", "D2", "D3"),
    "D3h": ("C2v", "C3", "C2", "C3h", "Cs", "C3v", "D2", "D3"),
    "Cinfv": ("K", "Cs"),
    "D7h": ("C2v", "C2", "C7", "C7v", "Cs", "C7h", "D7", "D2"),
    "Dinfh": ("C2v", "C2h", "C2", "K", "Kh", "Cs", "Cinfv", "Ci"),
    "S4": ("C2",),
    "C5h": ("C5", "Cs"),
    "D7d": ("C2v", "D2h", "C2h", "C2", "C7", "C7v", "Cs", "Ci", "D7", "D2"),
    "Th": (
        "C2v",
        "S6",
        "D2h",
        "C2h",
        "C3",
        "C2",
        "Cs",
        "T",
        "D3d",
        "C3v",
        "Ci",
        "D2",
        "D3",
    ),
    "Ci": (),
    "C7h": ("C7", "Cs"),
    "Ih": (
        "C2v",
        "S6",
        "D2h",
        "C2h",
        "C3",
        "C2",
        "C5",
        "C5v",
        "I",
        "D5d",
        "Cs",
        "T",
        "D3d",
        "C3v",
        "Th",
        "Ci",
        "D5",
        "D2",
        "D3",
    ),
    "D6": ("C3", "C2", "C6", "D2", "D3"),
    "D8": (),
    "Cs": (),
    "D7": ("C2", "C7", "D2"),
    "D4": ("C2", "C4", "D2"),
    "D5": ("C2", "C5", "D2"),
    "D2": ("C2",),
    "D3": ("C3", "C2", "D2"),
    "C1": (),
}


def init_hashing(depth, width, alg):
    """Initialize the static variables of the hashing module.

    The module is ManipulateAggregates.collection.hashIO. Please make sure that
    depth * width is not larger than the number of hexadecimal places of
    the returned hash values by the chosen algorithm alg.

    Args:
        depth: (int) number of subdirectories to be created
        width: (int) number of letters per subdirectory name
        alg: (string) hashing algorithm. See the function
            ManipulateAggregates.collection.hashIO.set_hashalg for
            supported algorithms.
    """
    hashIO.set_depth(depth)
    hashIO.set_width(width)
    hashIO.set_hashalg(alg)


def double_array(mylist):
    """Create a C array of doubles from a list.
    
    Args:
        mylist: (list of numerical values) will be converted to a C array
    """
    c = doubleArray(len(mylist))
    for i, v in enumerate(mylist):
        c[i] = v
    return c


def double_dist(iterable):
    """Returns a list of (array,-array) where array is a
    C-array made from a single element of iterable.
    
    Args:
        iterable: numpy array of (float,float,float) a numpy array of
            3D-vectors that are to be converted to plain C-arrays.
    """
    return [(double_array(i), double_array(-i)) for i in iterable]


def print_dx_file(prefix, hash, dictionary, values, comment):
    """Helper function to make writing a DX-file easier and to reduce the
    number of arguments that has to be passed from function to function.
    Uses hashing to reduce the number of files per directory.

    Args:
        prefix: (string, can contain a path) a prefix to the name of the
            to-be-printed dx file
        hash: (bool) whether or not to use hashing to reduce the number of
            files per directory
        dictionary: (dictionary with entries: filename, counts, org, delx,
            dely, delz, gzipped) describes the specifics of the to-be-printed
            dx file
        values: (list of floats) volumetric data that is stored in the dx file
        comment: (string) this will be the first line of the dx file but
            prefixed with a hash ('#')
    """
    if hash:
        filename = hashIO.hashpath(prefix + dictionary["filename"])
        directory = os.sep.join(filename.split(os.sep)[0:-1])
        if os.path.exists(directory):
            if not os.path.isdir(directory):
                raise FileExistsError(
                    errno.ENOTDIR,
                    "Cannot create necessary directory, file exists but no directory",
                    directory,
                )
        else:
            try:
                os.makedirs(directory)
            # catch a race condition occuring when several workers happen to create the same directory alost simultaneously
            except OSError as e:
                if e.errno == errno.EEXIST:
                    if not os.path.isdir(directory):
                        raise FileExistsError(
                            errno.ENOTDIR,
                            "Race condition likely, file exists but no directory",
                            directory,
                        )
                    else:
                        print(
                            "WARNING: Race condition likely, file exists and (luckily) is directory: %s"
                            % (directory),
                            file=sys.stderr,
                        )
                else:
                    raise e
    else:
        filename = prefix + dictionary["filename"]
    counts = dictionary["counts"]
    org = dictionary["org"]
    delx = dictionary["delx"]
    dely = dictionary["dely"]
    delz = dictionary["delz"]
    gzipped = dictionary["gzipped"]
    c_print_dx_file(
        filename,
        counts,
        org,
        delx,
        dely,
        delz,
        values,
        comment=comment,
        gzipped=gzipped,
        formatstring="14.13e",
    )


def no_none_string(string):
    """Determine whether a string is the literal 'None'.

    Args:
        string: (string) will be tested

    Returns:
        whether ot not the string is 'None'
    """
    return not (string == "None")


def general_grid(org, countspos, countsneg, dist, postprocessfunc=None, resetval=False):
    """Return a 3D-grid.

    Args:
        org: (NumPy array of shape (3,) and dtype float) the origin of the grid
        countspos: (NumPy array of shape (3,) and dtype int) how many points
            shall be taken in positive x,y and z directions
        countsneg: (NumPy array of shape (3,) and dtype int) same as
            countspos but for negative directions
        dist: (NumPy array of shape (3,) and dtype int) the distances for x,y
            and z directions
        postprocessfunc: (function) if not None, this function is applied to the
            grid prior to returning it
        resetval: (bool) if Ture, return (grid,reset_val) instead of just grid.
            Here, reset_val is the last element of the original grid that has been
            removed from the grid an can be used to return the vector to org.
            Meaning adding all elements in grid to org leaves the vector at org if
            resetval==False and at org-reset_val if resetval==True.

    Returns:
        the grnerated grid or (grid,reset_val), depending on the value of
        resetval.
    """
    # create helper vectors
    start = org - countsneg * dist
    end = org + countspos * dist
    # create grid for rotation of the molecule
    space = [
        numpy.linspace(s, e, num=cp + cn + 1)
        for s, e, cp, cn in zip(start, end, countspos, countsneg)
    ]
    # space     = [numpy.linspace(s,e,num=cp+cn+1,dtype=float)
    #             for s,e,cp,cn
    #             in zip(start,end,countspos,countsneg)
    #            ]
    a1 = numpy.array(
        [(x,) for x in space[0] for y in space[1] for z in space[2]], dtype=float
    )
    a2 = numpy.array(
        [(y,) for x in space[0] for y in space[1] for z in space[2]], dtype=float
    )
    a3 = numpy.array(
        [(z,) for x in space[0] for y in space[1] for z in space[2]], dtype=float
    )
    # a1,a2,a3  = numpy.array(numpy.meshgrid(*space,indexing="ij"))
    # a1.shape  = (-1,1)
    # a2.shape  = (-1,1)
    # a3.shape  = (-1,1)
    grid = numpy.concatenate((a1, a2, a3), axis=1)
    if postprocessfunc is not None:
        grid = postprocessfunc(grid)
    if resetval:
        reset = -grid[len(grid) - 1]
        return grid, reset
    else:
        return grid


def prepare_molecules(mol1, mol2, aligned_suffix="", save_aligned=False, align=True):
    """Prepare two molecules for the scanning procedure. First, align mol1
    and mol2 with their centers to org and their third and second main
    axes with [1,0,0] and [0,1,0], respectively.  Then, append mol2 to mol1 and
    enable tagging to allow easily moving mol2 with one command.  The
    OBAggregate object is returned.

    Args:
        mol1: (object of type ManipulateAggregates.aggregate.agg) the molecule
            that is to be kept fixed during the scan
        mol2: (object of type ManipulateAggregates.aggregate.agg) the molecule
            that is to be moved around during the scan
        aligned_suffix: (string) add this to the original file names to get
            the files to which the aligned geometries shall be saved
        save_aligned: (bool) whether or not to save the aligned geometries
        align: (bool) whether or not to perform the alignment procedure

    Returns:
        an OBAggregate object that contains mol1 and mol2
    """
    nr_scan_mols = mol2.obmol.GetNrMolecules()
    nr_fix_mols = mol1.obmol.GetNrMolecules()
    if align:
        # align the molecule's longest axis with the x-axis and the second longest axis
        # with the y-direction and center the molecule to the origin
        mol1.align([0, 0, 0], [1, 0, 0], [0, 1, 0])
        mol2.align([0, 0, 0], [1, 0, 0], [0, 1, 0])
        if save_aligned:
            mol1.write(mol1.info["name"] + aligned_suffix)
            mol2.write(mol2.info["name"] + aligned_suffix)
    else:
        mol1.obmol.Center()
        mol2.obmol.Center()
    # append the molecule
    mol1.append(mol2)
    obmol = mol1.obmol
    # configure the aggregate for using tags in order to be able to move multiple
    # molcules in the aggregate with one command
    obmol.EnableTags()
    tag = obmol.CreateTag(nr_scan_mols)
    for i in range(nr_fix_mols, nr_fix_mols + nr_scan_mols):
        obmol.AddToTag(i, tag)
    return obmol


def get_old_dxfiles(olddirs, suffix):
    """Search directories for previously created dx files.

    This function uses the module ManipulateAggregates.collection.hashIO to
    search in subdirectories with hashed names.

    Args:
        olddirs: (list of strings) the pathnames of the directories in which
            the previously generated dx-files are present.
        suffix: (string) the names of the dx files are I_SUFFIX.dx where I is
            an integer number and SUFFIX is the value of this parameter.

    Returns:
        a list of strings containing the pathnames of the old dx files.
    """
    olddxfiles = {}
    dxregex = re.compile("^[1-9][0-9]*_%s$" % (suffix))
    for d in set(olddirs):
        # ignore empty directory names
        if len(d) == 0:
            continue
        oldlength = len(olddxfiles)
        if os.path.isdir(d):
            for f in hashIO.listfiles(d, dxregex, nullsize=False, nulldepth=False):
                olddxfiles[int((f.split(os.sep)[-1]).split("_")[0])] = f
            if len(olddxfiles) == oldlength:
                print(
                    "WARNING: directory supposed to contain dx files from previous runs %s does not contain anything matching ^[1-9][0-9]*_%s$ . Skipping."
                    % (d, suffix),
                    file=sys.stderr,
                )
    return olddxfiles
