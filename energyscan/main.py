import sys
from ConfigParser import NoOptionError

import numpy as np

from collection.read import read_config_file as rf

class WrongJobtypeError(Exception):
    pass

def _print_example():
    """
    Print an example config file for an energyscan to stdout.
    """
    s="""#all lines starting with # are comments and can be removed

###GENERAL VALUES###
#declare the jobtype. NO DEFAULT SO IT MUST BE PROVIDED. Values are: scan, minimasearch
#May be a comma-separated list of jobtypes which will then be performed in sequence
jobtype         = scan,minimasearch
#use a cubic spatial grid that is not truncated. optional, default: full
sp_gridtype     = full
#declare how many steps in the positive x,y and z directions shall be used
countsxyz       = 50,50,50
#declare the stepsize in x,y and z directions
distxyz         = 0.5,0.5,0.5
#Whether or not to get progress reports during the calculation
#0: suppress progress reports
#1: print detailed progress reports during computation (SCAN: only works if OMP_NUM_THREADS is 1 or not set)
#2: print general progress reports (SCAN: whenever an angle was scanned, MINIMASEARCH: whenever minima for one angle were determined)
progress        = 2
#If True, only perform checks for the given config file but do not perform any computations. optional, default: False
config_check    = False
#IMPORTANT NOTICE: declare the exact same grid for a MINIMASEARCH jobtype that was used for a previous SCAN run!

###JOBTYPE SCAN###
geometry1       = aligned.xyz
#use the same geometry again (general variable replacement)
geometry2       = %(geometry1)s
#declare the force field. Select one of: mmff94, ghemical, uff, gaff. optional, default: mmff94
forcefield      = uff
#use a cubic angular grid that is not truncated. optional, default: full
ang_gridtype    = full
#declare how many steps in the positive angular directions of the main axes shall be used
countspos       = 1,1,1
#declare how many steps in the negative angular directions of the main axes shall be used
countsneg       = 0,0,0
#declare the stepsize in the directions of the main axes
dist            = 30.0,30.0,30.0
#if vdW surfaces are farther apart than this, do not evaluate the energy. optional, default: 100.0
cutoff          = 100.0
#ignore if <0.0 but if >0.0, scale all vdw-radii by this value before trying to determine clashes. optional, default: -1.0
vdw_scale       = -1.0
#True if dx-files shall be saved. optional, default: True
save_dx         = True
#how many columns shall be used in the dx file. optional, default: 3
columns         = 3
#the name of the dx-files (will be prepended by number and underscore). optional, default: out.dx
suffix          = out.dx
#do you want the aligned structures to be saved? optional, default: True
save_aligned    = True
#given the input geometryes, give the suffix for the aligned structures. optional, default: .aligned
aligned_suffix  = .aligned
#prefix this to any minimum energy geometry that will be saved. optional, default: template_
prefix          = template_
#do you want to save the global energy minima per angular arrangement? optional, default: True
save_noopt      = True
#do you want to save the global energy minima per angular arrangement aftger performing a force field optimization? optional, default: False
save_opt        = False
#steps for that force field optimization
optsteps        = 500
#This value must be larger than any other energy value
#you expect since all filtered values will be set to this
maxval          = 1000000000
#If True, after scanning all anergies, set all values that are
#at least 'maxval' to the true total maximum. optional, default: False
correct         = False
#if True, will find the optimum angle for every spatial arrangement. optinal, default: False
#all values in sp_opt files are in the exact same order
sp_opt          = True
#save the optimum energies to the following dx-file. optional, default: sp_opt.dx
#will be skipped if value None is given
sp_opt_dx       = sp_opt.dx
#save the corresponding geometries to the following xyz-file. optional, default: sp_opt.xyz
#will be skipped if value None is given
sp_opt_xyz      = sp_opt.xyz
#save the corresponding angles to the following csv file. optional, default: sp_opt_ang.csv
#will be skipped if value None is given
sp_opt_ang      = sp_opt_ang.csv
#like correct, but for the spatial grid. optional, default: True
sp_correct      = True
#do you want to remove such entries from sp_opt_ang and sp_opt_xyz where vdW clashes
#occured or where the molecules' vdW surfaces were farther apart from each other than curoff?
#If False, entries in the csv file that would be removed are given angles 360,360,360
#and entries in the xyz file will show two completely overlapping molecules. optional, default: False
sp_remove       = True
#do you want the global optimum to be saved to globalopt.xyz? optional, default: True
globalopt       = True

###JOBTYPE MINIMASEARCH###
#how to check whether two points on the grid are neighbours.
#Possible values (without quotes): 'eukledian' (Eukledian metric),'manhattan_single' (Manhattan metric with 
# one cutoff for all directions),'manhattan_multiple' (Manhattan metric with one cutoff for each direction)
#optional, default: manhattan_multiple
neighbour_check_type = manhattan_multiple
#distance below which two points are considered neighbours. Must be 'float,float,float' for 'manhattan_multiple'
# and float otherwise. optional, default: distxyz (same as one grid spacing in each direction)
distance_cutoff = %(distxyz)s
#a value by which the distance_cutoff will be scaled. That way, getting 2 shells of neighbours is as
# easy as setting this value to something slightly greater than 2. optional, default: 1.1
cutoff_scale    = 1.1
#Whether to compute and how to sort the depth-values of the minima. optional, default: 1
#Depth values are unused as of now, so you can as well declare a value of 0
#0: switch off sorting and do not determine depths of minima
#1: depths is equal to the difference of the value at the minimum and the average of all neighbours
#2: depths is equal to the difference of the value at the minimum and the minimum of all neighbours
depths_sort     = 1
#the minimum value that the volumetric data has to be lower than that of all its neighbours to be
# considered a minimum (may also be negative). optional, default: 0.0
degeneration    = 0.0
#how many neighbours do you want to search per gridpoint. optional, default: auto (works only
# for 'manhattan_multiple')
nr_neighbours   = auto
#how many neighbours do you expect a single point to have (at most). optional, default: nr_neighbours
# Greatly impacts performance, must be greater or equal nr_neighbours
max_nr_neighbours = %(nr_neighbours)s
#from where to take the volumetric data. optional, default: 'from_scan,.' Possible values:
#from_scan: take those dx-files that a jobtype of type scan would create if it had this config file (with adjusted jobtype).
#           I.e. 'from_scan,DIR' would take all dx-files created by a scan in the directory DIR. "." matches the current directory.
#dir_regex: take those dx-files that match the given regular expression. They will be sorted by the first integer number
#           in the name. I.e. 'dir_regex,/home/test/dir,\\\\.dx$' would match everything ending on ".dx" in "/home/test/dir".
#           Please double backslashes. The regular expression and DIR must not contain commas.
volumetric_data   = from_scan,.

###JOBTYPE SIMILARITYSCREENING###
#WARNING: NOT YET IMPLEMENTED, ONLY A DUMMY SECTION
#input data:
#   the file to which the information about the minima was saved
#       format: 1st line: nr of minima, 2nd line: comment line
#       3rd to 5th line: dx-file-like header for spatial grid
#       6th to 8th line: dx-file-like header for angular grid
#       then on each line: spatial_grid_index    x y z (coords)    angular_grid_index   a1 a2 a3 (angles)    value    depth (maybe)
#   how many geometries the user wants
#   maybe upper and lower bounds for RMSD and energy cutoffs, definitely steppsize for bineary search
"""
    print s

def _main(input_file):
    config = {
            "config_check"   : "False",
            "forcefield"     : "mmff94",
            "geometry1"      : "%(geometry)s",
            "geometry2"      : "%(geometry)s",
            "sp_gridtype"    : "full",
            "cutoff"         : "100.0",
            "vdw_scale"      : "-1.0",
            "ang_gridtype"   : "full",
            "save_dx"        : "True",
            "columns"        : "3",
            "suffix"         : "out.dx",
            "save_aligned"   : "True",
            "prefix"         : "template_",
            "aligned_suffix" : ".aligned",
            "save_noopt"     : "True",
            "save_opt"       : "False",
            "optsteps"       : "500",
            "progress"       : "2",
            "correct"        : "False",
            "maxval"         : "1000000000",
            "sp_opt"         : "False",
            "sp_opt_dx"      : "sp_opt.dx",
            "sp_opt_xyz"     : "sp_opt.xyz",
            "sp_opt_ang"     : "sp_opt_ang.csv",
            "sp_correct"     : "True",
            "sp_remove"      : "True",
            "globalopt"      : "True",
            "distance_cutoff": "%(distxyz)s",
            "cutoff_scale"   : "1.1",
            "degeneration"   : "0.0",
            "depths_sort"    : "1",
            "nr_neighbours"  : "auto",
            "volumetric_data": "from_scan,.",
            "neighbour_check_type" : "manhattan_multiple",
            "max_nr_neighbours"    : "%(nr_neighbours)s"
            }
    parser = rf(input_file,defaults=config)
    if parser.get_boolean("config_check"):
        print "This is a check of the config file."
    try:
        jobtype_list = parser.get_str("jobtype")
    except NoOptionError:
        raise KeyError("Necessary option 'jobtype' missing from config file")
    for jobtype in jobtype_list.split(","):
        if jobtype == "scan":
            from energyscan.scan import scan_main as jobtype_main
        elif jobtype == "minimasearch":
            from energyscan.minimasearch import minimasearch_main as jobtype_main
        else:
            raise WrongJobtypeError("Wrong jobtype chosen. Supported ones are: scan, minimasearch")
        print "Running jobtype %s..."%(jobtype)
        jobtype_main(parser)
        print "...finished jobtype %s"%(jobtype)
    if parser.get_boolean("config_check"):
        print "Config file seems fine."

if __name__ == "__main__":
    if len(sys.argv)==1:
        _print_example()
    else:
        for infile in sys.argv[1:]:
            _main(infile)
