"""This is the executable for the ManipulateAggregates.energyscan submodule.

Usage is as "energyscan [OPTIONS] CONFIGFILE1 [CONFIGFILE2] [...]"

Parallelization is supported. See the documentation for
ManipulateAggregates.energyscan for further details.
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
import sys
import operator

from ManipulateAggregates import get_data_dir
from ManipulateAggregates.collection.read import read_config_file
from ManipulateAggregates.collection.read import NoOptionInConfigFileError
from ManipulateAggregates.energyscan.scan import scan_main
from ManipulateAggregates.energyscan.minimasearch import minimasearch_main
from ManipulateAggregates.energyscan.similarityscreening import similarityscreening_main

# Default process name
PROCNAME = "EScan"
try:
    from FireDeamon import set_procname

    set_procname(PROCNAME)
except ImportError:
    set_procname = lambda s: None


class WrongJobtypeError(Exception):
    pass


global SHORTHELPTEXT
# The short help text message
SHORTHELPTEXT = r"""Usage:
    energyscan [OPTIONS] CONFIGFILE1 [CONFIGFILE2] [...] 

Command line OPTIONS:
    --help                  print this message
    --longhelp              print a long help message (which is also the default config
                            file). Comments in this message explain the meanings 
    --porphin-example       run an example scan for porphin in your current directory as
                            publised in the paper "A Program for Automatically
                            Predicting Supramolecular Aggregates and Its Application to
                            Urea and Porphin" by Sachse et al, accessible at
                            https://dx.doi.org/10.1002/jcc.25151 (beware: takes a very
                            long time and uses up more than 30GB of disk space) 
    --urea-example          run the urea scan from the same paper (beware: takes a long
                            time and uses up approx. 30GB of disk space) 
    --anthracene-example    run a quick scan using the anthracene molecule 
"""

global LONGHELPTEXT
## the long help text message (also a default config file)
LONGHELPTEXT = r"""# This is an example config file that also tries to give some explanations about what
# all the parameters do.
# Use the program as "energyscan CONFIGFILE1 [CONFIGFILE2] [...]
# Keywords are case-insensitive.
# All lines starting with # are comments and can be removed.

# ##VALUES NEEDED BY SEVERAL JOBTYPES AND GENERAL VALUES###
# declare the jobtype. NO DEFAULT SO IT MUST BE PROVIDED. Values are: scan, minimasearch
# May be a comma-separated list of jobtypes which will then be performed in sequence
jobtype         = scan,minimasearch,similarityscreening
# declare the gridtype to be used. Possible values are:
#    full: a cubic spatial grid that is not truncated
#    half: a cubic spatial grid that can be truncated to only comprise half-spaces
# optional, default: full
sp_gridtype     = full
# declare how many steps in the x,y and z directions shall be used for "full" and "half"
# grids. Note that for gridtype "half", the actual extent of the grid is modified by the
# value for "halfspace"
countsxyz       = 50,50,50
# declare the stepsize in x,y and z directions for "full" and "half" grids
distxyz         = 0.5,0.5,0.5
# declare whether the grid (only works for gridtypes "full" and "half" so far) shall
# automatically be adjusted to better fit the system at hand. Possible values are
# "countsxyz", "distxyz", "none" and empty depending on wether countsxyz or distxyz
# shall be adjusted or whether no adjustment shall be performed.  optional, default:
# EMPTY
sp_autoadjust   =
# declare the name of a file to which information about the grid shall be saved.
# optional, default: spgrid.dat
sp_gridsave     = spgrid.dat
# for gridtype "half", declare which half-spaces shall be used in x, y and z directions
# in the form of a vector of integers. A positive (negative) value indicates that the
# positive (negative) half-space shall be used. A value of 0 indicates that both halfs
# shall be used.  Hence, 0,0,0 is identical to gridtype "full". Additionaly, to avoid
# border effects, a value of e.g. n or -n causes n additional points to be considered in
# the neglected half-space (this is handy when using nr_neighbours = auto) optional,
# default: 0,0,0
halfspace       = 0,0,0
# Whether or not to get progress reports during the calculation
# 0: suppress progress reports
# 1: print detailed progress reports during computation (SCAN: only works if
#    OMP_NUM_THREADS is 1 or not set)
# 2: print general progress reports (SCAN: whenever an angle was scanned, MINIMASEARCH:
#    whenever minima for one angle were determined)
progress        = 2
# If True, only perform checks for the given config file but do not perform any
# computations. optional, default: False
config_check    = False
# which geometries to use for SCAN and SIMILARITYSCREENING. The default for geometry2 is
# to use the same geometry again (general variable replacement).
geometry1       =
geometry2       = %(geometry1)s
# whether or not you want to align the molecules with their center to 0,0,0 and their
# third/second main axis to the x/y-axis prior to any calculation. optional, default:
# True
prealign        = True
# whether or not you want to have dxfiles written in gzipped format to save disk space.
# This will put more load on the processor when reading or writing. optional, default:
# False
gzipped         = False
# IMPORTANT NOTICE: declare the exact same grid for a MINIMASEARCH jobtype that was used
# for a previous SCAN run!

# ##JOBTYPE SCAN###
# give I1/I2 (with I1 and I2 being integers and I2>=I1). If not equal 1/1, this job is
# part of a split job and will only perform a certain subset (i.e., partition) of the
# scan (e.g. 1/5 would perform the first fifth, 2/5 the second fifth). An exception will
# be raised for any jobtype other than scan if this value is not 1/1.  optional,
# default: 1/1
partition       = 1/1
# declare the force field. Select one of: mmff94, ghemical, uff, gaff. optional,
# default: mmff94
forcefield      = uff
# use a cubic angular grid that is not truncated. optional, default: full
ang_gridtype    = full
# declare how many steps in the positive angular directions of the main axes shall be
# used
countspos       = 1,1,1
# declare how many steps in the negative angular directions of the main axes shall be
# used
countsneg       = 0,0,0
# declare the stepsize in the directions of the main axes
dist            = 30.0,30.0,30.0
# if vdW surfaces are farther apart than this, do not evaluate the energy. optional,
# default: 100.0
cutoff          = 100.0
# ignore if <0.0 but if >0.0, scale all vdw-radii by this value before trying to
# determine clashes. optional, default: -1.0
vdw_scale       = -1.0
# True if dx-files shall be saved. optional, default: True
save_dx         = True
# how many columns shall be used in the dx file. optional, default: 3
columns         = 3
# the name of the dx-files (will be prepended by number and underscore). optional,
# default: out.dx
suffix          = out.dx
# do you want the aligned structures to be saved? optional, default: True
save_aligned    = True
# given the input geometryes, give the suffix for the aligned structures. optional,
# default: .aligned
aligned_suffix  = .aligned
# prefix this to any minimum energy geometry that will be saved. optional, default:
# template_
prefix          = template_
# do you want to save the global energy minima per angular arrangement? optional,
# default: True
save_noopt      = True
# do you want to save the global energy minima per angular arrangement after performing
# a force field optimization? optional, default: False
save_opt        = False
# steps for that force field optimization
optsteps        = 500
# This value must be larger than any other energy value you expect (in units of the
# selected forcefield) since all filtered values will be set to this
maxval          = 1000000000
# If True, after scanning all energies, set all values that are at least 'maxval' to the
# true total maximum. optional, default: False
correct         = False
# if True, will find the optimum angle for every spatial arrangement. optinal, default:
# False all values in sp_opt files are in the exact same order
sp_opt          = True
# save the optimum energies to the following dx-file. optional, default: sp_opt.dx will
# be skipped if value None is given
sp_opt_dx       = sp_opt.dx
# save the corresponding geometries to the following xyz-file. optional, default:
# sp_opt.xyz will be skipped if value None is given
sp_opt_xyz      = sp_opt.xyz
# save the corresponding angles to the following csv file. optional, default:
# sp_opt_ang.csv will be skipped if value None is given
sp_opt_ang      = sp_opt_ang.csv
# like correct, but for the spatial grid. optional, default: True
sp_correct      = True
# do you want to remove such entries from sp_opt_ang and sp_opt_xyz where vdW clashes
# occured or where the molecules' vdW surfaces were farther apart from each other than
# curoff?  If False, entries in the csv file that would be removed are given angles
# 360,360,360 and entries in the xyz file will show two completely overlapping
# molecules. optional, default: False
sp_remove       = True
# do you want the global optimum to be saved to globalopt.xyz? optional, default: True
globalopt       = True
# if this is a restarted scan, declare the directories (comma separated) where all
# previous data can be found. optional, default: EMPTY
scan_restartdirs =
# to limit the number of files per directory, this program uses the first letters of the
# hash of a dx-file's name to put it in subdirectories. Declare the hashing algorithm to
# be used (get a list of all supported ones via python -c 'import hashlib;print
# hashlib.algorithms;'), string, optional ,default: md5
hashalg         = md5
# declare how many hex-digits shall be used per directory level, optional, default: 2
hashwidth       = 2
# declare how many levels of directories shall be used for the hasing process, optional,
# default: 2
hashdepth       = 2

# ##JOBTYPE MINIMASEARCH###
# how to check whether two points on the grid are neighbours.
# Possible values (without quotes): 'eukledian' (Eukledian metric),'manhattan_single'
# (Manhattan metric with one cutoff for all directions),'manhattan_multiple' (Manhattan
# metric with one cutoff for each direction) optional, default: manhattan_multiple
neighbour_check_type = manhattan_multiple
# distance below which two points are considered neighbours. Must be 'float,float,float'
# for 'manhattan_multiple' and float otherwise. The special value 'auto', which also
# considers auto-adjustments of the spatial grid, is only supported for
# 'manhattan_multiple' and grids 'full' or 'half'. optional, default: 'auto'
distance_cutoff = auto
# a value by which the distance_cutoff will be scaled. That way, getting 2 shells of
# neighbours is as easy as setting this value to something slightly greater than 2.
# optional, default: 1.1
cutoff_scale    = 1.1
# Whether to compute and how to sort the depth-values of the minima. optional, default:
# 1 Depth values are unused as of now, so you can as well declare a value of 0 0: switch
# off sorting and do not determine depths of minima 1: depths is equal to the difference
# of the value at the minimum and the average of all neighbours 2: depths is equal to
# the difference of the value at the minimum and the minimum of all neighbours
depths_sort     = 1
# the minimum value that the volumetric data has to be lower than that of all its
# neighbours to be considered a minimum (may also be negative). optional, default: 0.0
degeneration    = 0.0
# how many neighbours do you want to search per gridpoint. optional, default: auto
# (works only for 'manhattan_multiple')
nr_neighbours   = auto
# how many neighbours do you expect a single point to have (at most). optional, default:
# nr_neighbours Greatly impacts performance, must be greater or equal nr_neighbours
max_nr_neighbours = %(nr_neighbours)s
# from where to take the volumetric data. optional, default: 'from_scan,.' Possible
# values:
#
# from_scan: take those dx-files that a jobtype of type scan would create if it had this
#     config file (with adjusted jobtype).  I.e. 'from_scan,DIR' would take all dx-files
#     created by a scan in the directory DIR. "." matches the current directory.  This
#     option respects the value of 'scan_restartdirs' and will also use those dx-files.
#
# dir_regex: take those dx-files that match the given regular expression. They will be
#     sorted by the first integer number in the name. I.e.
#     'dir_regex,/home/test/dir,\\.dx$' would match everything ending on ".dx" in
#     "/home/test/dir".  Please double backslashes. The regular expression and DIR must
#     not contain commas.
volumetric_data   = from_scan,.
# declare the file to which the data about the minima shall be saved. If ending in '.gz'
# (without quotes), it will be be gzipped. optional, default: minima.dat
minima_file_save  = minima.dat
# have the worker processes wait for the main processes after they finished their jobs
# for the declared number of angular arrangements. This allows the main process to keep
# up in cases where the worker processes are too fast. Try decreasing this if the main
# process shows too-high memory usage. optional, default: 100
pool_chunksize    = 100

# ##JOBTYPE SIMILARITYSCREENING###
# from where to take the data about the minima that were found. If ending in '.gz'
# (without quotes), it will be considered to be gzipped. optional, default: same as
# minima_file_save
minima_file_load  = %(minima_file_save)s
# how many geometries the user wants at least. Those geometries are as diverse as
# possible in their geometries.  Will try to find the number closest to the given one,
# but you might also get fewer depending on the cutoffs for RMSD and energy.
nr_geometries     = 10
# the maximum RMSD-cutoff for SIMILARITYSCREENING.
rmsd_min          = 1
# starting from rmsd_max, increase the RMSD-cutoff by this value until fewer than
# nr_geometries were found.  Suppose that one is called rmsd_max, return the geometries
# for rmsd_max minus rmsd_step.
rmsd_step         = 0.5
# only consider geometries whose energy is closer to that of the global minimum geometry
# for the screenig.  A negative value switches off screening by energy. This value is
# given in meV (milli electron volts).  optional, default: -100
energy_cutoff     = -100
# whether or not to use the given energy cutoff in force field units (differs with force
# field). optional, default: False
use_ff_units      = False
# declare the xyz-file to which all geometries shall be saved that passed the screening
# by RMSD and possibly energy. optional, default screened.xyz
screened_xyz      = screened.xyz
# declare how many decimal places shall be used when screening for symmetry equivalence.
# A value smaller than 0 turns off this behaviour. optional, default: 2 (highly
# recommended, decrease if you get duplicates)
symprec           = 2
# whether or not to align the obtained conformers after the simmilarity screening.
# optional, default: True
postalign         = True
# declare the maximum number of screening steps that are to be performed (to aviod
# infinite loops).  optional, default: 500
maxscreensteps    = 500
# whether or not to determine (and possibly screen by) pointgroups of the conformers.
# optional, default: False
pointgroups       = False
# whether or not to include subgroups. If True, a conformer with a higher symmetry is
# also considered to be part of all pointgroups that are subgroups of the determined
# pointgroups (e.g., C4 conformers would also belong to the C2 pointgroup). optional,
# default: False
subgroups         = False
# whether or not to always exclude the C1 pointgroup. optional, default: True
exclude_c1        = True
# after which similarity screening step to perform the determination of pointgroups. The
# special keywords "first" and "last" are accepted. Otherwise, an integer >=0 must be
# provided. optional, default: first
pgstep            = first
# if True, save all conformers to separate files depending on their pointgroups (see
# "pgfile" suffixed with the pointgroup name). If False, perform only screening by
# pointgroup (based on pgregex). If pgregex is not given and this is False, the
# pointgroup determination is skipped. optional, default: True
pgwrite           = True
# declare the prefix for the xyz files to which to save the conformers belonging to the
# pointgroups. If the given value already ends on a known extension, this will be
# removed prior to suffixing.
pgfile            = %(geometry1)s
# declare a Python regular expression. All pointgroups whose name is not matched by this
# regular expression are screened. Beware: you must double backslashes! If, for
# instance, you want to exclude everything but C2, use "^(?!C2$)" (without quotes). If
# you want to exclude all pointgroups starting with D, use "^(?!D.*$)" optional,
# default: EMPTY (matches everything, i.e., no screening, BEWARE the value of
# exclude_c1)
pgregex           =
# declare a comma-separated list of indices (counting starts at 1) that indicate
# hydrogens that shall not be discarded during the similarity screening (i.e.,
# computation of RMSD, determination of equivalence). optional, default: EMPTY (i.e.,
# there are no important hydrogens)
consider_h1       =
# the same as the above for the second molecule. If the special value "SAME" is given
# and geometry1 is equal to geometry2, the same indices will be used for the second
# molecule. An error will be raised if geometry1 is unequal geometry2, consider_h1 is
# not empty and consider_h2 is "SAME".
consider_h2       = SAME
"""


def _print_example():
    """Print an example config file for an energyscan to stdout."""
    print(LONGHELPTEXT)


def _print_help():
    """Print the help message"""
    print(SHORTHELPTEXT)


global DEFAULT_CONFIG
## default config options
DEFAULT_CONFIG = {
    "aligned_suffix": ".aligned",
    "ang_gridtype": "full",
    "columns": "3",
    "config_check": "False",
    "consider_h1": "",
    "consider_h2": "SAME",
    "correct": "False",
    "cutoff": "100.0",
    "cutoff_scale": "1.1",
    "degeneration": "0.0",
    "depths_sort": "1",
    "distance_cutoff": "auto",
    "energy_cutoff": "-100",
    "exclude_c1": "True",
    "forcefield": "mmff94",
    "geometry2": "%(geometry1)s",
    "globalopt": "True",
    "gzipped": "False",
    "halfspace": "0,0,0",
    "hashalg": "md5",
    "hashdepth": "2",
    "hashwidth": "2",
    "max_nr_neighbours": "%(nr_neighbours)s",
    "maxscreensteps": "500",
    "maxval": "1000000000",
    "minima_file_load": "%(minima_file_save)s",
    "minima_file_save": "minima.dat",
    "neighbour_check_type": "manhattan_multiple",
    "nr_neighbours": "auto",
    "optsteps": "500",
    "package_data_dir": get_data_dir(),
    "partition": "1/1",
    "pgfile": "%(geometry1)s",
    "pgregex": "",
    "pgstep": "first",
    "pgwrite": "True",
    "pointgroups": "False",
    "pool_chunksize": "100",
    "postalign": "True",
    "prealign": "True",
    "prefix": "template_",
    "progress": "2",
    "save_aligned": "True",
    "save_dx": "True",
    "save_noopt": "True",
    "save_opt": "False",
    "scan_restartdirs": "",
    "screened_xyz": "screened.xyz",
    "sp_autoadjust": "",
    "sp_correct": "True",
    "sp_gridsave": "spgrid.dat",
    "sp_gridtype": "full",
    "sp_opt_ang": "sp_opt_ang.csv",
    "sp_opt_dx": "sp_opt.dx",
    "sp_opt": "False",
    "sp_opt_xyz": "sp_opt.xyz",
    "sp_remove": "True",
    "subgroups": "False",
    "suffix": "out.dx",
    "symprec": "2",
    "use_ff_units": "False",
    "vdw_scale": "-1.0",
    "volumetric_data": "from_scan,.",
}

global MANDATORY_OPTIONS
# Mandatory options for certain jobtypes.
#
# The following options have to be provided in the config file for the
# following jobtypes:
#
#  - jobtype: always
#  - countsxyz: scan, minimasearch
#  - distxyz: scan, minimasearch
#  - geometry: scan, similarityscreening
#  - countspos: scan, minimasearch
#  - countsneg: scan, minimasearch
#  - dist: scan
#  - nr_geometries: similarityscreening
#  - rmsd_min: similarityscreening
#  - rmsd_step: similarityscreening
MANDATORY_OPTIONS = {
    "scan": [
        "countsneg",
        "countspos",
        "countsxyz",
        "dist",
        "distxyz",
        "geometry1",
        "geometry2",
        "jobtype",
    ],
    "minima search": ["countsneg", "countspos", "countsxyz", "distxyz", "jobtype",],
    "similarity screening": [
        "geometry1",
        "geometry2",
        "jobtype",
        "nr_geometries",
        "rmsd_min",
        "rmsd_step",
    ],
}


def _main(input_file):
    # default config
    config = DEFAULT_CONFIG
    options = [o for o in config] + list(
        set([mo for mopts in MANDATORY_OPTIONS.values() for mo in mopts])
    )
    parser = read_config_file(input_file, defaults=config, nocase=True)
    jobtype_list = parser.get_str("jobtype")
    # jobtypes have long and short names but both shall be treated the same so the
    # following is a mapping of the long and short forms to a unified form
    jobtype_dict = {
        "scan": "scan",
        "s": "scan",
        "minimasearch": "minima search",
        "ms": "minima search",
        "similarityscreening": "similarity screening",
        "ss": "similarity screening",
    }
    functions_dict = {
        "scan": scan_main,
        "minima search": minimasearch_main,
        "similarity screening": similarityscreening_main,
    }
    # Jobs have to be performed in a certain order to make sense. The order is given
    # in the following dictionary (starting at 0 and increasing):
    order_dict = {"scan": 0, "minima search": 1, "similarity screening": 2}
    try:
        # sort the jobtypes to be performed by an increasing order parameter
        jobtype_list = sorted(
            # since each job must not be performed multiple times, use a set to get the
            # irreducible set of jobtypes
            set(
                # create a list of tuples each consisting of the desired jobtype name
                # and its order parameter
                [
                    (jobtype_dict[e.lower()], order_dict[jobtype_dict[e.lower()]])
                    for e in jobtype_list.split(",")
                ]
            ),
            key=operator.itemgetter(1),
        )
    except KeyError as e:
        raise ValueError(
            "Given short or long form does not match any known jobtype: %s" % (e)
        )
    # check whether all mandatory options are present
    missing_options = []
    for jobtype, discard in jobtype_list:
        for opt in MANDATORY_OPTIONS[jobtype]:
            try:
                parser.get_str(opt)
            except NoOptionInConfigFileError:
                missing_options.append(opt)
    if len(missing_options) > 0:
        print(
            "ERROR: could not find the following mandatory options in the config file:",
            file=sys.stderr,
        )
        for o in missing_options:
            print(o, file=sys.stderr)
        raise NoOptionInConfigFileError("Incomplete input.")
    del missing_options
    unknown_options = parser.check_against(options)
    if len(unknown_options) > 0:
        print(
            "WARNING: the following are unknown lines in the config file:",
            file=sys.stderr,
        )
        for o in unknown_options:
            print(o, file=sys.stderr)
        print(file=sys.stderr)
    del unknown_options
    if parser.get_boolean("config_check"):
        print("This is a check of the config file.")
    for jobtype, __ in jobtype_list:
        jobtype_main = functions_dict[jobtype]
        print("Running %s..." % (jobtype))
        jobtype_main(parser)
        set_procname(PROCNAME)
        print("...finished %s\n" % (jobtype))
    if parser.get_boolean("config_check"):
        print("Config file seems fine.")


def entrypoint():
    if len(sys.argv) == 1:
        _print_example()
        print("\n\nHelp message:")
        _print_help()
    else:
        for arg in sys.argv[1:]:
            if arg == "--help":
                _print_help()
            elif arg == "--longhelp":
                _print_example()
            elif arg in ("--porphin-example", "--urea-example", "--anthracene-example"):
                molecule = arg.split("-")[2]
                sample_script = os.path.join(get_data_dir(), "{}.cfg".format(molecule))
                geometry = os.path.join(get_data_dir(), "{}.xyz".format(molecule))
                print(
                    "Running example scan using config file: {}".format(sample_script)
                )
                print("Using geometry file: {}".format(geometry))
                _main(sample_script)
            else:
                _main(arg)


if __name__ == "__main__":
    entrypoint()
