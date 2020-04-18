# -*- coding: utf-8 -*-
"""This is the executable for manipulating aggregates.
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
import io
import copy
import re

import ManipulateAggregates as ma
from maagbel import pybel

global use_np
try:
    import numpy

    use_np = True
except ImportError:
    use_np = False

# default process name
PROCNAME = "ManipAgg"
try:
    from FireDeamon import set_procname

    set_procname(PROCNAME)
except ImportError:
    set_procname = lambda s: None

global FUNCTIONDICT, VSDICT, CPDICT, AGGREGATES, CURRENTAGG, FF, CONFORMER, FORMAT, ENVIRONMENTS, SET, TAGGING, PART
AGGREGATES = []
CURRENTAGG = None
ENVIRONMENTS = []
SET = True
TAGGING = False

global HELPTEXT, RENDERHELPTEXT, POTHELPTEXT, VISHELPTEXT, MANIPHELPTEXT, AUXHELPTEXT
global CLOSERHELPTEXT
# the short help text message
HELPTEXT = """This script manipulates internal degrees of freedom of a molecule or aggregate.

The custom reduced version of OpenBabel (https://github.com/razziel89/MaAg-bel) is
required. Supports all filetypes supported by that reduced version (you can run the
command `manipagg --list formats` to get a list). The default input filetype is guessed
from the extension. The output filetype is the input filetype by default. Much
functionality requires libFireDeamon (https://github.com/razziel89/libFireDeamon).

Usage (switches are positional, i.e., affect everything behind them):
    manipagg [GEOMETRYFILE] [SWITCHES]

Some switches require a GEOMETRYFILE, some others don't.

Mandatory options for long versions are mandatory for short versions, too. You can
separate switch and value by either space or the equals sign. A "#" following a swich
means it requires a number of parameters separated by spaces. E.g., --switch #2 means
that this switch requires 2 arguments which must be separated by any number of spaces >0.
The symbol "[#X]" means that X parameters are optional but if the first is given, all the
others have to be given, too.

Please note that parameters starting with dashes are fine unless they are the same as an
option. In that case, you cannot use it.

Command line switches:
~~~~~~~~~~~~~~~~~~~~~~

--help|-h       Print this help and exit
--render-help   See more detailed help about automatically rendering a visualization
--pot-help      See more detailed help about how to define an electrostatic potential
                and how to obtain it
--vis-help      See more detailed help about how to modify visualizations
--manip-help    See more detailed help about how to manipulate a geometry
--aux-help      See information about some auxilliary switches.
--full-help     See all help texts (useful for grep'ing for switches)
--ff #1         Declare force field to use ("None" is also possible, switching
                off everything that requires one) (default: mmff94)
--intype|-i #1  Set the type of the input file (default: guess from filename)
--infile|-I #1  Give the name of the input file (not required if the input file is the
                first argument)
--outtype|-o #1 Set the type of the output file (default: guess from filename)
--outfile|-O #1 Set the name of the output file (default: do not output anything)
--conf #1       Declare which conformer from a file that can contain
                multiple conformers (such as the xyz-format) you wish to load
--list [#1]     List supported plugin options. To get a list of plugins, pass 'plugins'
                as argument to this switch or pass no argument
--example-vdw   Run an example visualization of the electrostatic potential on a
                molecule's van-der-Waals surface as publised in the paper "Introducing
                double polar heads to highly fluorescent Thiazoles: Influence on
                supramolecular structures and photonic properties" by Kaufmann et al,
                accessible at https://doi.org/10.1016/j.jcis.2018.04.105
--example-iso   Run an example visualization of the same molecule as --example-vdw on an
                isosurface of the molecule's electron density. This will put some files
                in your current directory.


"""
## help message for molecule manipulation
MANIPHELPTEXT = """More information about geometry manipulation switches:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

--bond|-b #1        Set the length of the bond defined by 2 atom indices to the desired
                    value in Anstroms, e.g., set the bond between atoms 1 and to to 5.5
                    Angstroms: 1,2=5.5 . If one atom is marked with a star (i.e. a
                    preceeding `*` without the tics) it will be kept fixed. Otherwise,
                    both atoms wil be moved halfway.
--angle|-a #1       Set an angle defined by 3 atom indices to a desired value in degrees,
                    e.g., 1,2,3=90.5 .
--dihedral|-d       Set the dihedral angle defined by 4 atom indices to the desired value
                    in degrees. A cis configuration corresponds to 0 degrees, a trans
                    configuration to 180 degrees, e.g., 1,2,3,4=90
--get|-g            Instead of setting the desired internal degrees of freedom the script
                    outputs them (ignores all angle or bondlength values given, e.g.
                    --bond 1,2=5 would result in the bondlength defined by atoms 1 and 2
                    to be output and the 5 will be ignored.)
--set|-s            Unset --get for everything following --set
--write             Do an intermediate write out of the output file to the specified file.
                    Very handy if things should be done in succession
--tag               UNDOCUMENTED
--app #1            Append the given second molecule to the current one. Uses the current
                    input file type. If options are given after --app and before --end,
                    they will be applied to the to-be-appended molecule.
--dup               Like --app, but it does not take an argument. Instead, the
                    current molecule is being duplicated. This switch also
                    requires a matching --end.
--end|--licate      Declare that no more actions are to applied to the to be appended
                    molecule. Both can be used with --app and --dup.
--gl #1             Glue the given second molecule to the one given. Behaves just as
                    --app with respect to --ue instead of --end.
--ue #2             Declare that no more actions are to applied to the to be glued
                    molecule.  As arguments, give two pairs of indices as i1,i2 and
                    m1,m2. Those pairs declare pairs of atoms the bonds between which
                    are to be cut. Then, the molecules will be glued together and i1 and
                    m1 will remain in the molecule and i2 and m2 (and all atoms connected
                    to them) will be cleaved.
--translate|-t #1   Translate the molecule by the given vector, e.g., --translate=1,0,0
                    would translate the molecule by 1 Angstrom in the x-direction
--rotate|-r #1      Rotate the molecule arount the given axis by the given angle e.g.
                    --rotate=1,0,0=90 would rotate the molecule by 90 degrees around the
                    x-axis
--rotate-main #1    Rotate the molecule around the given main axis, e.g.,
                    --rotate-main=1=90 would rotate the molecule by 90 degrees around its
                    first main axis
--mirror #2         Mirror the molecule at a plane or a point (inversion). Declare the
                    normal vector first and a point in the plane next. If the former is
                    0,0,0, point inversion will be performed. Example: 1,0,0 0,0,0 would
                    mirror at the yz-plane and 0,0,0 0,0,0 would invert at the origin.
--mirror-center #2  The same as --mirror, but perform all these actions after centering
                    to 0,0,0. It will be moved back to original center afterwards.
--align #3          Align the molecule with its center to a given point and align the
                    third and second main axes to the two axes given e.g. --align 0,0,0
                    1,0,0 0,1,0 would align the molecule's center to the origin and its
                    third/second main axis with the x/y-axis The third main axis is
                    usually the longest extent.
--part [#1]         Apply all subsequent manipulateions (translation, rotation, alignment,
                    mirroring) only to the specified covalently bound subunit (usually
                    a molecule). Counting starts at 0. Leave out options to switch back to
                    treating everyting together.
--cleave #1         Cleave one part of a molecule by cutting a bond in half. That bond is
                    specified as an argument consisting of two indices separated by a
                    comma.
--optimize #1       Perform a force-field optimization with the given number of steps
                    (default: 500)
--closer #2 [#1]    Move two parts of an aggregate closer together with respect to their
                    centers until a vdW-clash occurs. Give as ---move-closer p1,p2
                    s (f,a) See --closer-help for additional information.
--closer-vec #3 [#1]
                    Move two parts of an aggregate closer together in the direction of
                    the vector given until a vdW-clash occurs or the distance between the
                    centers increases. Give as ---closer-vec p1,p2 v1,v2,v3 s (f,a) See
                    --closer-help for additional information.

"""
## help message for molecule visualization
VISHELPTEXT = """More information about visualization switches:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

--visualize-pot|--vp #1 [#1] [#1]
                    Visualize the molecule and the potential on a van-der-Waals surface
                    using OpenGL, Python, libFireDeamon and CGAL. It will automatically
                    be aligned with its longest axis from right to left.
                    You specify: Z [R] [F1,F2]
                      - Z is a zoom factor (float)
                      - R is the number of refinement steps (int).
                      - F1 and F2 (floats) are the scaling factor and the shrink factor
                            for the skin-surface generation, respectively.
--visualize-iso|--vi #2 [#1] [#1]
                    Like --visualize-pot, but plot the electrostatic potential on an isosurface
                    instead of a vdW surface.
                    You specify: Z I [L] [C]
                      - Z is a zoom factor (float)
                      - I is the iso value (default: 0.005)
                      - L is a comma-separated list of atom indices (ints) around which
                            the isosurface shall be built. Special values: all, noH, auto
                            (auto is the default)
                      - C is P1,P2,A,R:
                        - P1 is the first of CGAL's surface mesh precisions (float)
                        - P2 is the second of CGAL's surface mesh precisions (float)
                        - A is the minimum angle between surface facets (float)
                        - R is the relative surface mesh generation precision (float)
--visualize-simple|--vs #1 [#1]
                    Visualize the molecule as spheres with corresponding
                    vdW-radii using OpenGL and Python.
                    You specify: Z S
                      - Z is a zoom factor (float)
                      - S is a scaling factor for the vdW spheres (float)
--window-title|--title #1
                    Set the title of the visualization window, which is also the prefix
                    for images saved to disk.
--window-resolution|--resolution|--res #1
                    Set the resolution of the visualization window as x,y (two ints).
--hide              Do not show the OpenGL window (useful for rendering a renderpath)
--swap-align        Usually, the molecule's third main axis is aligned perpendicular to
                    the visualization plane, its second main axis is aligned to the
                    left-right direction and it's center of mass is moved to the center
                    of the screen. This command suspends and re-enables this alignment
                    procedure (first occurence disables, second occurence enables, etc.)
--contrast #1       By supplying "high" as a parameter (default), the color scale is:
                    blue (negative) - black (vanishing) - red (positive) for the
                    visualization of the electrostatic potential.
                    By supplying "low" as a parameter, the color scale is: red (negative)
                    - yellow (less negative) - green (vanishing) - turquoise (less
                    positive) - blue (positive) for the visualization of the
                    electrostatic potential.
--invert            Invert potential data no matter where it has been obtained from
--svgscale #1       Save an SVG file with the given name that shows the color scale.
--save-vis #2       When visualizing, save the visualization state. You specify: W F
                      - W is an arbitrary combination of the words start and end
                            stating whether you want the visualization state saved at the
                            beginning or the end, respectively. "None" turns that off.
                      - F is the name of the file to which to save the sate.
                            A prefix might be added.
                    Press comma during visualization to save additional visualization
                    states. Does not work for --visualize-simple.
--load-vis #1       Load visualization data from the given file. Will also initiate
                    visualization.
--povray #1         Declare an integer. If this integer is >0, the resolution of an image
                    rendered using PoVRay will have this times the resolution of the
                    OpenGL window. If the integer is <=0, support for PoVRay is switched
                    off (default: 1, <=0 not recommended).
--povlight #1 [#1]
                    Declare an axis (three comma-separated floats) and an angle
                    (in degrees) that define a rotation for all normal vectors
                    prior to PoVRay visualization. If the axis is "frontal",
                    "front" or "straight", illumination will happen directly
                    from the front. The default is to slightly rotate all
                    normal vectors as this looks nicer.
--refscale #2 [#n]  You provide: R D1 [D2] [...]
                      - R is a Python regular expression that will be used to match
                            against files in the given directories
                      - D1 is a directory (as are all other DN)
                    The color scale of the potential plot will be adjusted so that all
                    scales, defined in the save files whose names match the regular
                    expression in the given directory, fall within the same overall scale
                    (to make them comparable). Incompatible with --colorscale.
--colorscale #1 [#1]
                    You provide C1 [C2]:
                      - C1 is a special keyword (see below) or the float value used as
                            the lower end of the color scale
                      - C2 is the float value used as the upper end of the color scale
                            (ignored if C1 is a special keyword)
                    Special values are: "auto" (default), or "independent" (only first
                    letter checked), which causes the use of independent color scales for
                    positive and negative values.  The special value "dependent" (same as
                    independent) causes the use of the same color scale for positive and
                    negative values.

"""
## auxilliary help message
AUXHELPTEXT = """Auxilliary switches that only output information:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

--dipole-moment #1  Output the molecule's dipole moment as obtained from the
                    specified charge method. See "manipagg --list charges" for a list of
                    available methods
--energy            Output the energy of the molecule according to the current force
                    field
--rmsd #1           Output RMSD between the current molecule and the given one as
                    well as the maximum difference in a single coordinate and for an
                    atom. This will use the currently defined intype and perform an
                    alignment beforehand to exclude influences from translation or
                    rotation.
--vdw-check         Output "True" if no two atoms from different molecules are closer
                    together than the sum of their vdW-radii and "False" otherwise
--spinmultiplicity  Output the molecule's spinmultiplicity
--pbond|--pb #2     Output the length of the bond defined by 2 atom indices in Angstroms
                    projected onto a given vector (2nd atgument).
--hlb               quick estimation of a molecule's HLB value
--repickle          if a visualization state cannot be loaded, that migh tbe because the
                    state was saved in Python 2 and you try to load it using Python 3.
                    This will try to convert the state to a more compatible representation.
                    WARNING: the original file will be overwritten!
--pgroup|--pg       Print the point group of the given structure.

"""
## help message for moving molecules closer together
CLOSERHELPTEXT = """Help text for the --closer command and the --closer-vec command:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WARNING: if one geometry file that was read in contains multiple geometries, that has to
be considered!

--closer #2 [#1]    Move two parts of an aggregate closer together with respect to their
                    centers until a vdW-clash occurs. Give as ---move-closer p1,p2
                    s (f,a)
--closer-vec #3 [#2]
                    Move two parts of an aggregate closer together in the direction of
                    the vector given until a vdW-clash occurs or the distance between the
                    centers increases. Give as ---closer-vec p1,p2 v1,v2,v3 s (f,a)

p1 and p2: 
            indices, starting at 0, indicating the molecules in the aggregate that shall
            be moved closer together.
v1,v2,v3: 
            components of the vector in which the first molecule shall be moved (will be
            inverted for the second one). 
s:          stepsize for movement (good value: 0.2). 
f:          factor by which all vdW-radii will be multiplied (default: 0.9). 
a:          value that is added to all vdW-radii (default: 0.0). 

"""
## help message for rendering a visualization path
RENDERHELPTEXT = """Help text for the --renderpath command:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The --renderpath command It takes one argument.

A renderpath tells the script how to manipulate the visualization of the given molecule.
All images that are rendered can be saved to disk.

A simple trajectory string looks as follows (spaces only to emphasize logical groups):

chain_of_commands | chain_of_values_separated_by_dashes / number_of_frames

The above can be repeated as often as desired if separated by commas. Apart from the
first command, each chain of commands has to follow a comma. You have to declare as many
values (separated by dashes) as you have declared commands.

Commands that can be chained:
    r1+:    rotate positively around first axis 
    r1-:    rotate negatively around first axis 
    r2+:    rotate positively around second axis 
    r2-:    rotate negatively around second axis 
    r3+:    rotate positively around third axis 
    r3-:    rotate negatively around third axis 
    t1+:    translate positively along first axis 
    t1-:    translate negatively along first axis 
    t2+:    translate positively along second axis 
    t2-:    translate negatively along second axis 
    t3+:    translate positively along third axis 
    t3-:    translate negatively along third axis 
    z+:     increase zoom level (default zoom level: 10) 
    z-:     decrease zoom level (default zoom level: 10) 

Special commands that do not take values or a number of frames and have to be the last
ones in the trajectory:
    n:      Do not save OpenGL images to disk 
    p:      Render every image via PoVRay 
    d:      Drop to an interactive view first where the user can rotate
            the molecule using the keybord. After a press of ESC, the
            path will be followed 
    s:      At each image, save the visualization state to disk. Requires --save-vis to
            be set 

Values: 
    rotation:    angles in degrees 
    translation: lengths in the unit given in the geometry file (usually Angstroms) 
    zoom:        change in zoom level 

Number of frames: 
The number of frames during which the given change in visualization will
be performed. A linear mapping from change to frame number is applied.

Example: 
`r1+r2-t3+z-|180-90-2-5/100,t1-z+|1-2/200,n,d`

First, the user will see the molecule and have the opportunity to change the view by
using the keyboard. After pressing ESC, the trajectory will start: In the first 100
frames, rotate around the first axis by 180° whilst rotating negatively around the second
axis by 90°, translating the molecule by 2 along the third axis and reducing the zoom
level by 5. In the next 200 frames, traslate negatively along the first axis by 1 whilst
increasing the zoom level by 2. None of the frames will be rendered to disk.

"""
## help message for obtaining an electrostatic potential or charges
POTHELPTEXT = """More information about switches regarding potentials, charges and densities:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please note that --density only works together with --orbitals.

--charges #1    Compute the potential from existing volumetric charge data
--empirical #1  Compute the electrostatic potential from charges obtained via an
                empirical method. You specify the method (see "obabel -L charges").
                Default is "mmff94".
--inter #1 [#3] Compute the potential by interpolating existing data.
                You specify: M [C1 C2 C3]
                  - M is the interpolation method ("distance" or "nearest" for inverse
                        distance weighting or nearest neighbour interpolation (default))
                  - C1 is the expotential parameter
                  - C2 is the root
                  - C3 is the cutoff (negative value switches this off)
                C2, C2 and C3 are only used for inverse distance weighting.
--orbitals      Compute the density or potential from molecular orbital data
--absolute      Charges at the atomic sites are absolute ones (opposite of --partial)
--partial       Charges at the atomic sites are partial ones (default)
--cube #1       Specify a CUBE file as the input file for volumetric data
--cube-vis #1   Specify a CUBE file as the input file for volumetric data (for isosurface
                generation, overwrites a previous --cube for this purpose)
--dx #1         Specify a DX file as the input file for volumetric data
--dx-vis #1     Specify a DX file as the input file for volumetric data (for isosurface
                generation, overwrites a previous --dx for this purpose)
--molden #1     Specify a Molden file as the input file for orbital data
--xyz #1        Specify a XYZ file as the input file for volumetric data. This type
                of file has Cartesian coordinates in the first 3 columns followed by the
                value at that point in space.
--density       The property to compute is the electron density. Only useful with --grid.
--potential     The property to compute is the electrostatic potential (default).
--grid [#1] [#1]
                Compute the specified property on a grid. You specify: P F
                  - P is the number of points in each of the 3 Cartesian directions for
                      the regular grid (default: 100)
                  - F is the file of type DX to which the data shall be saved (default:
                      potential.dx or density.dx depending on the property)

"""
## help message for default key bindings for visualization
KEYSHELPTEXT = """Key bindings for the visualization window:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Key bindings: 
    ESC : quit 
    = : zoom in 
    - : zoom out 
    w : move molecule up 
    s : move molecule down 
    a : move molecule left 
    d : move molecule right 
    q : move molecule to front 
    e : move molecule to back 
    i : rotate molecule positively around 1st axis 
    k : rotate molecule negatively around 1st axis 
    j : rotate molecule positively around 2nd axis 
    l : rotate molecule negatively around 2nd axis 
    u : rotate molecule positively around 3rd axis 
    o : rotate molecule negatively around 3rd axis 
    . : save OpenGL snapshot of current view 
    , : save current visualization for later restore 
    p : render approximation of current view via PoVRay 

"""


def _le(arg, l):
    if len(arg) != l:
        raise ValueError("Wrong format for argument.")


def _nn(arg):
    return arg is not None


def _vs():
    if CURRENTAGG is None:
        return VSDICT.__setitem__
    else:
        return CURRENTAGG.set_vs


def _gvs():
    if CURRENTAGG is None:
        return VSDICT.__getitem__
    else:
        return CURRENTAGG.get_vs


def _cp():
    if CURRENTAGG is None:
        return CPDICT.__setitem__
    else:
        return CURRENTAGG.set_cp


def _gcp():
    if CURRENTAGG is None:
        return CPDICT.__getitem__
    else:
        return CURRENTAGG.get_cp


def __help():
    print(HELPTEXT)
    sys.exit(0)


def __hlb():
    hlb, discard1, discard2 = ma.aggregate.hlb.compute(CURRENTAGG)
    print(hlb)


def __aux_help():
    print(AUXHELPTEXT)
    sys.exit(0)


def __closer_help():
    print(CLOSERHELPTEXT)
    sys.exit(0)


def __full_help():
    print(HELPTEXT)
    print(MANIPHELPTEXT)
    print(POTHELPTEXT)
    print(VISHELPTEXT)
    print(AUXHELPTEXT)
    print(RENDERHELPTEXT)
    print(CLOSERHELPTEXT)
    print(KEYSHELPTEXT)
    sys.exit(0)


def __manip_help():
    print(MANIPHELPTEXT)
    sys.exit(0)


def __pot_help():
    print(POTHELPTEXT)
    sys.exit(0)


def __render_help():
    print(RENDERHELPTEXT)
    sys.exit(0)


def __renderpath(path):
    _vs()("renderpath", path)


def __vis_help():
    print(VISHELPTEXT)
    print(KEYSHELPTEXT)
    sys.exit(0)


def __absolute():
    _cp()("partial", False)


def __add_h():
    CURRENTAGG.obmol.AddHydrogens()


def __align(center, axis3, axis2):  # 3
    center = list(map(float, center.split(",")))
    axis3 = list(map(float, axis3.split(",")))
    axis2 = list(map(float, axis2.split(",")))
    CURRENTAGG.align(center, axis3, axis2, part=PART)


def __angle(arg):  # 1
    if SET:
        temp = arg.split("=")
        _le(temp, 2)
        atoms, angle = temp
        angle = float(angle)
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 3)
        CURRENTAGG.set_angle(atoms[0], atoms[1], atoms[2], angle)
    else:
        temp = arg.split("=")
        atoms = temp[0]
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 3)
        print(CURRENTAGG.get_angle(atoms[0], atoms[1], atoms[2]))


def __app(filename):  # 1
    global ENVIRONMENTS, AGGREGATES, CURRENTAGG
    ENVIRONMENTS.append("append")
    AGGREGATES.append(CURRENTAGG)
    CURRENTAGG = None
    __infile(filename)


def __bond(arg):  # 1
    if SET:
        temp = arg.split("=")
        _le(temp, 2)
        atoms, angle = temp
        angle = float(angle)
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 2)
        CURRENTAGG.set_bondlength(atoms[0], atoms[1], angle)
    else:
        temp = arg.split("=")
        atoms = temp[0]
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 2)
        print(CURRENTAGG.get_bondlength(atoms[0], atoms[1]))


def __charges(arg):  # 1
    _cp()("type", "charges")


def __cleave(arg):  # 1
    atoms = list(map(int, arg.split(",")))
    _le(atoms, 2)
    CURRENTAGG.cleave(atoms[0], atoms[1])


def __closer(parts, stepsize, config=None):  # 2 [#1]
    parts = list(map(int, parts.split(",")))
    _le(parts, 2)
    step = float(stepsize)
    if _nn(config):
        config = list(map(float, config.split(",")))
        _le(config, 2)
        CURRENTAGG.move_closer(
            parts[0],
            parts[1],
            stepsize=step,
            vdw_factor=config[0],
            vdw_added=config[1],
            vec=None,
        )
    else:
        CURRENTAGG.move_closer(parts[0], parts[1], stepsize=step)


def __closer_vec(parts, vector, stepsize, config=None):  # 3 [#1]
    parts = list(map(int, parts.split(",")))
    _le(parts, 2)
    step = float(stepsize)
    vector = list(map(float, vector.split(",")))
    _le(vector, 3)
    if _nn(config):
        config = list(map(float, config.split(",")))
        _le(config, 2)
        CURRENTAGG.move_closer(
            parts[0],
            parts[1],
            stepsize=step,
            vdw_factor=config[0],
            vdw_added=config[1],
            vec=vector,
        )
    else:
        CURRENTAGG.move_closer(parts[0], parts[1], stepsize=step, vec=vector)


def __colorscale(start, end=None):  # 1 [#1]
    try:
        lowerscale = float(start)
    except ValueError:
        if start.startswith("d"):
            start = "dependent"
        elif start.startswith("i"):
            start = "independent"
        if start != "auto":
            _vs()("colorscale", start)
        return
    if _nn(end):
        upperscale = float(end)
        try:
            import cPickle as p
        except ImportError:
            import pickle as p
        scalefile = "SCALE_%f_%f" % (lowerscale, upperscale)
        d = {"face_colourscale": (lowerscale, upperscale)}
        f = io.open(scalefile, "wb")
        p.dump(d, f, 2)
        f.close()
        _vs()("colorscale", ("^%s$" % (scalefile), "."))
    else:
        raise TypeError(
            "The first argument is no special keyword but no second argument given."
        )


def __conf(conf):  # 1
    global CONFORMER
    if conf in ("first", "last"):
        CONFORMER = conf
    else:
        CONFORMER = int(conf)


def __contrast(c):  # 1
    if c.lower().startswith("h"):
        _vs()("high_contrast", True)
    elif c.lower().startswith("l"):
        _vs()("high_contrast", False)
    else:
        raise ValueError("--contrast only accepts 'high' or 'low'.")


def __cube(filename):  # 1
    _cp()("chargefiletype", "cube")
    _cp()("potfiletype", "cube")
    _vs()("isofiletype", "cube")
    _cp()("chargefiletype", "cube")
    _cp()("chargefile", filename)
    _cp()("potfile", filename)
    _vs()("isofile", filename)
    _cp()("chargefile", filename)


def __cube_vis(filename):  # 1
    _vs()("isofiletype", "cube")
    _vs()("isofile", filename)


def __density():
    _cp()("property", "density")


def __dihedral(arg):
    if SET:
        temp = arg.split("=")
        _le(temp, 2)
        atoms, angle = temp
        angle = float(angle)
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 4)
        CURRENTAGG.set_dihedral(atoms[0], atoms[1], atoms[2], atoms[3], angle)
    else:
        temp = arg.split("=")
        atoms = temp[0]
        atoms = list(map(int, atoms.split(",")))
        _le(atoms, 4)
        print(CURRENTAGG.get_dihedral(atoms[0], atoms[1], atoms[2], atoms[3]))


def __dipole_moment(method):  # 1
    oldmethod = CURRENTAGG.get_cp("method")
    _cp()("method", method)
    print(CURRENTAGG.get_dipole_moment())
    _cp()("method", oldmethod)


def __dup():  # 1
    global ENVIRONMENTS, AGGREGATES, CURRENTAGG
    ENVIRONMENTS.append("append")
    AGGREGATES.append(CURRENTAGG)
    CURRENTAGG = CURRENTAGG.duplicate()


def __dx(filename):  # 1
    _cp()("chargefiletype", "dx")
    _cp()("potfiletype", "dx")
    _vs()("isofiletype", "dx")
    _cp()("chargefiletype", "dx")
    _cp()("chargefile", filename)
    _cp()("potfile", filename)
    _vs()("isofile", filename)
    _cp()("chargefile", filename)


def __dx_vis(filename):  # 1
    _vs()("isofiletype", "dx")
    _vs()("isofile", filename)


def __empirical(method):  # 1
    _cp()("type", "empirical")
    _cp()("method", method)


def __end():
    global ENVIRONMENTS, AGGREGATES, CURRENTAGG
    if not ENVIRONMENTS.pop() == "append":
        raise ValueError("The switch --end has to follow --app.")
    newagg = CURRENTAGG
    CURRENTAGG = AGGREGATES.pop()
    CURRENTAGG.append(newagg)


def __energy():
    print(CURRENTAGG.get_energy())


def __ff(ff):  # 1
    global FF
    FF = ff


def __get():
    global SET
    SET = False


def __gl(filename):  # 1
    global ENVIRONMENTS, AGGREGATES, CURRENTAGG
    ENVIRONMENTS.append("glue")
    AGGREGATES.append(CURRENTAGG)
    CURRENTAGG = None
    __infile(filename)


def __grid(points="100", filename=None):  # [#1] [#1]
    if not use_np:
        raise RuntimeError("Could not import numpy, cannot perform a grid computation.")
    points = int(points)
    coordinates = numpy.array(CURRENTAGG.get_coordinates())
    min_corner = numpy.amin(coordinates, axis=0) - 10.0
    max_corner = numpy.amax(coordinates, axis=0) + 10.0
    counts_xyz = numpy.array([points, points, points])
    org_xyz = min_corner
    # grid creation copied from energyscan.scan but slightly altered
    space = [
        numpy.linspace(s, e, num=c)
        for s, e, c in zip(min_corner, max_corner, counts_xyz)
    ]
    # just take the difference between the first elements in every direction to get the stepsize
    delta_x = numpy.array([space[0][1] - space[0][0], 0.0, 0.0])
    delta_y = numpy.array([0.0, space[1][1] - space[1][0], 0.0])
    delta_z = numpy.array([0.0, 0.0, space[2][1] - space[2][0]])
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
    if _gcp()("property") == "density":
        if not _nn(filename):
            filename = "density.dx"
        prop = CURRENTAGG.get_density(grid)
    else:
        if not _nn(filename):
            filename = "potential.dx"
        prop = CURRENTAGG.get_potential(grid)
    ma.collection.write.print_dx_file(
        filename, counts_xyz, org_xyz, delta_x, delta_y, delta_z, prop
    )


def __hide():
    _vs()("hide", True)


def __infile(filename):  # 1
    global CURRENTAGG, VSDICT
    if CURRENTAGG is None:
        CURRENTAGG = ma.aggregate.read_from_file(
            filename, fileformat=FORMAT, conf_nr=CONFORMER, ff=FF
        )
        for k in VSDICT:
            CURRENTAGG.set_vs(k, VSDICT[k])
        VSDICT.clear()
        for k in CPDICT:
            CURRENTAGG.set_cp(k, CPDICT[k])
        CPDICT.clear()
    else:
        raise ValueError("You already specified a primary input file.")
    return


def __inter(method, exp=None, root=None, cutoff=None):  # 1 [#3]
    _cp()("type", "interpolation")
    _cp()("interpolation", method)
    if _nn(exp):
        _cp()("int_exponent", exp)
    if _nn(root):
        _cp()("int_root", root)
    if _nn(cutoff):
        _cp()("cutoff", cutoff)


def __intype(informat):  # 1
    global FORMAT
    FORMAT = informat


def __invert():
    _cp()("invert_potential", True)


def __list(type="plugins"):  # 1
    configs = pybel.getpluginconfigs(type)
    print("Configs for plugin '{}': {}".format(type, ", ".join(configs)))


def __load_vis(filename):  # 1
    _vs()("savestart", False)
    _vs()("saveend", False)
    ma.aggregate.visualize.RenderExtern(filename, agg=CURRENTAGG, dictionary=VSDICT),


def __mirror(point, normal):  # 2
    point = list(map(float, point.split(",")))
    normal = list(map(float, normal.split(",")))
    _le(point, 3)
    _le(normal, 3)
    CURRENTAGG.mirror(normal, point, center_it=False, part=PART)


def __mirror_center(point, normal):  # 2
    point = list(map(float, point.split(",")))
    normal = list(map(float, normal.split(",")))
    _le(point, 3)
    _le(normal, 3)
    CURRENTAGG.mirror(normal, point, center_it=True, part=PART)


def __molden(filename):  # 1
    _cp()("orbfiletype", "molden")
    _cp()("orbfile", filename)


def __optimize(steps):  # 1
    CURRENTAGG.optimize(int(steps))


def __orbitals():
    _cp()("type", "orbitals")


def __outfile(filename):  # 1
    CURRENTAGG.write(filename)


def __outtype(outtype):  # 1
    CURRENTAGG.info["outformat"] = outtype


def __part(part=None):
    global PART
    if part is not None:
        PART = int(part)
    else:
        PART = part


def __partial():
    _cp()("partial", True)


def __pbond(bond, vector):  # 2
    atoms = list(map(int, bond.split(",")))
    _le(atoms, 2)
    vector = list(map(float, vector.split(",")))
    _le(vector, 3)
    print(CURRENTAGG.get_bondlength(atoms[0], atoms[1], projection=vector))


def __pgroup():  # 0
    print(CURRENTAGG.get_pointgroup())


def __potential():
    _cp()("property", "potential")


def __povray(scale):  # 1
    _vs()("povray", int(scale))


def __povlight(axis, angle=None):  # 2
    if axis.lower() in ("frontal", "front", "straight"):
        axis, angle = ([1.0, 0.0, 0.0], 0.0)
    else:
        axis = axis.split(",")
        _le(axis, 3)
        axis = list(map(float, axis))
        angle = float(angle)
    _vs()("visrotmat", (axis, angle))


def __refscale(regex, dir1, *dirs):  # 2 [#n]
    dirs = [dir1] + list(dirs)
    _vs()("colorscale", (regex, "|".join(dirs)))


def __repickle(filename):  # 1
    try:
        import cPickle as p
    except ImportError:
        import pickle as p

        print("WARNING: cPickle module should be present in Python 2.", file=sys.stderr)
        print(
            "         Are you sure you are running this using Python 2?",
            file=sys.stderr,
        )
    f = io.open(filename, "rb")
    try:
        obj = p.load(f)
    except UnicodeDecodeError as e:
        print(
            "ERROR during unpickling. Maybe you did not use the same Python version as for pickling?",
            file=sys.stderr,
        )
        raise e
    f.close()
    for key in (
        "povray_data",
        "faces",
    ):
        if key not in obj:
            raise ValueError(
                "The loaded file %s is most likely no saved visualization state."
                % (filename)
            )
        try:
            obj[key] = [a.tolist() for a in obj[key]]
        except AttributeError:
            raise ValueError(
                "The loaded file %s has most likely already been repickled."
                % (filename)
            )
    f = io.open(filename, "wb")
    # protocol version 2 stays compatible with Python 2 (but is slower than more recent versions)
    p.dump(obj, f, 2)
    f.close()


def __rmsd(filename):  # 1
    global AGGREGATES, CURRENTAGG
    AGGREGATES.append(CURRENTAGG)
    CURRENTAGG = None
    __infile(filename)
    newagg = CURRENTAGG
    CURRENTAGG = AGGREGATES.pop()
    CURRENTAGG.rmsd(newagg, True)


def __rotate(arg):  # 1
    temp = arg.split("=")
    _le(temp, 2)
    axis, angle = temp
    angle = float(angle)
    axis = list(map(float, axis.split(",")))
    _le(axis, 3)
    CURRENTAGG.rotate(axis, angle, part=PART)


def __rotate_center(arg):  # 1
    if PART is not None:
        center = CURRENTAGG.obmol.GetCenterPart(PART)
        center = [-center.GetX(), -center.GetY(), -center.GetZ()]
    else:
        center = [-c for c in CURRENTAGG.get_center()]
    CURRENTAGG.translate(center, part=PART)
    __rotate(arg)
    CURRENTAGG.translate([-c for c in center], part=PART)


def __rotate_main(arg):  # 1
    temp = arg.split("=")
    _le(temp, 2)
    axis, angle = temp
    angle = float(angle)
    axis = int(axis)
    CURRENTAGG.rotate_main(axis, angle, part=PART)


def __rotate_main_center(arg):  # 1
    if PART is not None:
        center = CURRENTAGG.obmol.GetCenterPart(PART)
        center = [-center.GetX(), -center.GetY(), -center.GetZ()]
    else:
        center = [-c for c in CURRENTAGG.get_center()]
    CURRENTAGG.translate(center, part=PART)
    __rotate_main(arg)
    CURRENTAGG.translate([-c for c in center], part=PART)


def __save_vis(words, filename=None):  # 2
    if words.startswith("start") or words.endswith("start"):
        _vs()("savestart", True)
    if words.startswith("end") or words.endswith("end"):
        _vs()("saveend", True)
    if words.lower() in ("none", "never", "no"):
        _vs()("savestart", False)
        _vs()("saveend", False)
        if filename is not None:
            _vs()("savefile", filename)
    else:
        if filename is None:
            raise TypeError(
                "__save_vis() takes exactly 2 arguments when the first is not 'none'."
            )
        _vs()("savefile", filename)


def __svgscale(filename):
    _vs()("svgscale", filename)


def __set():
    global SET
    SET = True


def __spinmultiplicity():
    print(CURRENTAGG.obmol.GetTotalSpinMultiplicity())


def __swap_align():
    _vs()("align", not (_gvs()("align")))


def __tag(*args):
    global TAGGING
    if len(args) == 0:
        if TAGGING:
            CURRENTAGG.tag_parts(-1)
        else:
            CURRENTAGG.tag_parts(1)
        TAGGING = not (TAGGING)
    else:
        args = list(map(int, args))
        CURRENTAGG.tag_parts(args)
        TAGGING = True


def __translate(arg):  # 1
    arg = list(map(float, arg.split(",")))
    _le(arg, 3)
    CURRENTAGG.translate(arg, part=PART)


def __ue(pair1, pair2):  # 2
    global ENVIRONMENTS, AGGREGATES, CURRENTAGG
    pair1 = list(map(int, pair1.split(",")))
    pair2 = list(map(int, pair2.split(",")))
    if not ENVIRONMENTS.pop() == "glue":
        raise ValueError("The switch --ue has to follow --gl.")
    newagg = CURRENTAGG
    CURRENTAGG = AGGREGATES.pop()
    CURRENTAGG.glue(newagg, pair1[0], pair1[1], pair2[0], pair2[1])


def __vdw_check(scale=None):
    if _nn(scale):
        print(CURRENTAGG.vdw_check(scale))
    else:
        print(CURRENTAGG.vdw_check())


def __visualize_iso(zoom, iso=None, atoms=None, config=None):  # 2 [#1] [#1]
    _vs()("type", "iso")
    _vs()("zoom", float(zoom))
    if _nn(iso):
        _vs()("isovalue", float(iso))
    if _nn(atoms):
        temp = atoms.split(",")
        if len(temp) == 1:
            try:
                int(temp[0])
                _vs()("iso_atoms", int(temp[0]))
            except ValueError:
                _vs()("iso_atoms", temp[0])
        else:
            _vs()("iso_atoms", list(map(int, temp)))
    if _nn(config):
        temp = config.split(",")
        _le(temp, 4)
        _vs()("mesh_criteria", list(map(float, [temp[1], temp[2], temp[0]])))
        _vs()("rel_precision", float(temp[3]))
    CURRENTAGG.visualize()


def __visualize_pot(zoom, refine=None, factors=None):  # 1 [#1] [#1]
    _vs()("type", "vdw")
    _vs()("zoom", float(zoom))
    if _nn(refine):
        _vs()("refine", int(refine))
    if _nn(factors):
        temp = factors.split(",")
        _le(temp, 2)
        _vs()("vdw_scale", float(temp[0]))
        _vs()("shrink_factor", float(temp[1]))
    CURRENTAGG.visualize()


def __visualize_simple(zoom, scale=None):  # 1 [#1]
    _vs()("type", "simple")
    _vs()("zoom", float(zoom))
    if _nn(scale):
        _vs()("vdw_scale", float(scale))
    CURRENTAGG.visualize()


def __window_resolution(res):  # 1
    res = RESOLUTIONS.get(res.lower(), list(map(int, res.split(","))))
    _le(res, 2)
    _vs()("resolution", res)


def __window_title(title):  # 1
    _vs()("title", title)


def __write(filename):
    CURRENTAGG.write(filename)


def __xyz(filename):  # 1
    _cp()("chargefiletype", "xyz")
    _cp()("potfiletype", "xyz")
    _vs()("isofiletype", "xyz")
    _cp()("chargefiletype", "xyz")
    _cp()("chargefile", filename)
    _cp()("potfile", filename)
    _vs()("isofile", filename)
    _cp()("chargefile", filename)


def __example_iso():
    # Determine files to use for this example
    data_dir = ma.get_data_dir()
    molden_file = os.path.join(data_dir, "dye5.molden")
    print(
        r"""Running example visualization on an electron density iso-surface.
This will effectively execute the following command:

    manipagg -I $FILE --orbitals --density --molden $FILE --grid 100 dens.dx \
             --dx-vis dens.dx --potential --save-vis start vis_iso.masave \
             --visualize-iso 1.0 0.001

After this computation finishes, run the following command to visualize the already
computed visualization without recomputation:

    manipagg --load-vis start_vis_iso.masave

Remember that you can press "p" to render your current view via PoVRay if you opted to
install it.  Unfortunately, the OpenGL view cannot be translated perfectly to PoVRay,
which means you should zoom out a bit to avoid clipping the sides.

For $FILE, we use:
"""
    )
    print("    " + molden_file + "\n\nKeybindings follow.\n\n")
    print(KEYSHELPTEXT)
    # Execute commands in order
    __infile(molden_file)
    __orbitals()
    __density()
    __molden(molden_file)
    __grid(points=100, filename="dens.dx")
    __dx_vis("dens.dx")
    __potential()
    __save_vis("start", filename="vis_iso.masave")
    __visualize_iso(1.0, iso=0.001)


def __example_vdw():
    # Determine files to use for this example
    data_dir = ma.get_data_dir()
    molden_file = os.path.join(data_dir, "dye5.molden")
    print(
        r"""Running example visualization on a van-der-Waals surface.
This will effectively execute the following command:

    manipagg -I $FILE --orbitals --molden $FILE \
             --potential --save-vis start vis_vdw.masave \
             --visualize-pot 1.0

After this computation finishes, run the following command to visualize the already
computed visualization without recomputation:

    manipagg --load-vis start_vis_vdw.masave

Remember that you can press "p" to render your current view via PoVRay if you opted to
install it.  Unfortunately, the OpenGL view cannot be translated perfectly to PoVRay,
which means you should zoom out a bit to avoid clipping the sides.

For $FILE, we use:
"""
    )
    print("    " + molden_file + "\n\nKeybindings follow.\n\n")
    print(KEYSHELPTEXT)
    # Execute commands in order
    __infile(molden_file)
    __orbitals()
    __molden(molden_file)
    __potential()
    __save_vis("start", filename="vis_vdw.masave")
    __visualize_pot(1.0)


global RESOLUTIONS
## special resolution keywords accepted by this script and their resolutions
RESOLUTIONS = {
    "vga": (640, 480),
    "svga": (800, 600),
    "qhd": (960, 540),
    "wsvga": (1024, 600),
    "xga": (1024, 768),
    "xga+": (1152, 864),
    "wxga": (1280, 720),
    "wxga": (1280, 768),
    "wxga": (1280, 800),
    "sxga": (1280, 1024),
    "wxga": (1440, 900),
    "uxga": (1600, 1200),
    "wsxga+": (1680, 1050),
    "fhd": (1920, 1080),
    "hd": (1920, 1080),
    "wuxga": (1920, 1200),
    "wqhd": (2560, 1440),
    "wqxga": (2560, 1600),
}

## this dictionary associates each switch with a function that performs an operation
FUNCTIONDICT = {
    "--absolute": __absolute,
    "--add-h": __add_h,
    "-h": __add_h,
    "--align": __align,
    "--angle": __angle,
    "-a": __angle,
    "--app": __app,
    "--aux-help": __aux_help,
    "--bond": __bond,
    "-b": __bond,
    "--charges": __charges,
    "--cleave": __cleave,
    "--closer": __closer,
    "--closer-help": __closer_help,
    "--closer-vec": __closer_vec,
    "--colorscale": __colorscale,
    "--conf": __conf,
    "--contrast": __contrast,
    "--cube": __cube,
    "--cube-vis": __cube_vis,
    "--density": __density,
    "--dihedral": __dihedral,
    "-d": __dihedral,
    "--dipole-moment": __dipole_moment,
    "--dup": __dup,
    "--dx": __dx,
    "--dx-vis": __dx_vis,
    "--empirical": __empirical,
    "--end": __end,
    "--energy": __energy,
    "--example-vdw": __example_vdw,
    "--example-iso": __example_iso,
    "--ff": __ff,
    "--full-help": __full_help,
    "--get": __get,
    "-g": __get,
    "--gl": __gl,
    "--grid": __grid,
    "--help": __help,
    "-h": __help,
    "--hlb": __hlb,
    "--hide": __hide,
    "--infile": __infile,
    "-I": __infile,
    "--inter": __inter,
    "--intype": __intype,
    "-i": __intype,
    "--invert": __invert,
    "--licate": __end,
    "--list": __list,
    "--load-vis": __load_vis,
    "--manip-help": __manip_help,
    "--mirror": __mirror,
    "--mirror-center": __mirror_center,
    "--molden": __molden,
    "--optimize": __optimize,
    "--orbitals": __orbitals,
    "--outfile": __outfile,
    "-O": __outfile,
    "--outtype": __outtype,
    "-o": __outtype,
    "--part": __part,
    "--partial": __partial,
    "--pbond": __pbond,
    "--pb": __pbond,
    "--pgroup": __pgroup,
    "--pg": __pgroup,
    "--potential": __potential,
    "--pot-help": __pot_help,
    "--povray": __povray,
    "--povlight": __povlight,
    "--refscale": __refscale,
    "--render-help": __render_help,
    "--renderpath": __renderpath,
    "--repickle": __repickle,
    "--rmsd": __rmsd,
    "--rotate-main": __rotate_main,
    "--rotate-main-center": __rotate_main_center,
    "--rotate": __rotate,
    "--rotate-center": __rotate_center,
    "-r": __rotate,
    "--save-vis": __save_vis,
    "--svgscale": __svgscale,
    "--set": __set,
    "-s": __set,
    "--spinmultiplicity": __spinmultiplicity,
    "--swap-align": __swap_align,
    "--tag": __tag,
    "--translate": __translate,
    "-t": __translate,
    "--ue": __ue,
    "--vdw-check": __vdw_check,
    "--vis-help": __vis_help,
    "--visualize-iso": __visualize_iso,
    "--vi": __visualize_iso,
    "--visualize-pot": __visualize_pot,
    "--vp": __visualize_pot,
    "--visualize-simple": __visualize_simple,
    "--vs": __visualize_simple,
    "--window-resolution": __window_resolution,
    "--resolution": __window_resolution,
    "--res": __window_resolution,
    "--window-title": __window_title,
    "--title": __window_title,
    "--write": __write,
    "--xyz": __xyz,
}

# the default conformer to be read in
CONFORMER = 1
# the default force field to use
FF = "mmff94"
# the default file format (None means: guess)
FORMAT = None
# store all visualization options until an input molecule was specified
VSDICT = {}
# store all porential options until an input molecule was specified
CPDICT = {}
# which covalently bound entity shall be treated
PART = None


def _expand_tilde(filename):
    if filename.startswith("~"):
        homedir = filename.split(os.sep)[0]
        if homedir == "~":
            filename = re.sub("^~", os.environ["HOME"], filename)
        else:
            username = homedir.split("~")[1]
            homedir = os.sep.join(os.environ["HOME"].split(os.sep)[:-1])
            filename = re.sub("^~", homedir + os.sep, filename)
    return filename


def _parse_commandline(argv):
    while len(argv) > 0:
        temp = argv.pop().split("=")
        switch = temp[0]
        if switch.startswith("--"):
            switch = switch.lower()
        try:
            func = FUNCTIONDICT[switch]
        except KeyError as e:
            raise ValueError("Switch %s not known, aborting." % (switch))
        if len(temp) == 1:
            sargs = ""
            fargs = []
        else:
            sargs = "=".join(temp[1:])
            fargs = ["=".join(temp[1:])]
        while len(argv) > 0 and not argv[-1].split("=")[0] in FUNCTIONDICT:
            temp = argv.pop()
            sargs += " %s" % (temp)
            fargs.append(temp)
        fargs = [_expand_tilde(f) for f in fargs]
        yield (switch, sargs, func, fargs)


def entrypoint():
    argv = sys.argv
    argv.reverse()
    argv.pop()
    if len(argv) == 0:
        raise ValueError("No arguments provided.")
    # treat special case of geometry file being the first argument
    if not argv[-1] in FUNCTIONDICT and os.path.isfile(argv[-1]):
        argv.append("--infile")
    for switch, sargs, func, fargs in _parse_commandline(argv):
        try:
            func(*fargs)
        except TypeError as e:
            print(
                "You probably specified a wrong number of arguments for the switch %s. Arguments: %s. Stacktrace follows."
                % (switch, sargs),
                file=sys.stderr,
            )
            raise e
        except AttributeError as e:
            print(
                "You probably did not specify a primary input file, but switch %s requires one. Stacktrace follows."
                % (switch),
                file=sys.stderr,
            )
            raise e
        except IndexError as e:
            print(
                "You probably specified --end or --ue without specifying --app or --gl beforehand. Your switch was: %s. Stacktrace follows."
                % (switch),
                file=sys.stderr,
            )
            raise e


if __name__ == "__main__":
    entrypoint()
