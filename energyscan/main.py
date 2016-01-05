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
jobtype         = scan
#declare the jobtype. NO DEFAULT SO MUST BE PROVIDED. Values are: scan
#May be a comma-separated list of jobtypes
geometry1       = aligned.xyz
#use the same geometry again (general variable replacement)
geometry2       = %(geometry1)s
#declare the force field. Select one of: mmff94, ghemical, uff, gaff. optional, default: mmff94
forcefield      = uff
#use a cubic grid that is not truncated. optional, default: full
sp_gridtype     = full
#declare how many steps in the positive x,y and z directions shall be used
countsxyz       = 50,50,50
#declare the stepsize in x,y and z directions
distxyz         = 0.5,0.5,0.5
#use a cubic grid that is not truncated. optional, default: full
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
#given the input geometryes, give the suffic for the aligned structures. optiona, default: .aligned
aligned_suffix  = .aligned
#prefix this to any minimum energy geometry that will be saved. optional, default: template_
prefix          = template_
#do you want to save the global energy minima per angular arrangement? optional, default: True
save_noopt      = True
#do you want to save the global energy minima per angular arrangement aftger performing a force field optimization? optional, default: False
save_opt        = False
#steps for that force field optimization
optsteps        = 500
#0: suppress progress reports
#1: print progress reports during computation (only works if OMP_NUM_THREADS is 1 or not set)
#2: print progress reports whenever an angle was scanned
progress        = 2
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
globalopt       = True"""
    print s

def _main(input_file):
    config = {
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
            "globalopt"      : "True"
            }
    parser = rf(input_file,defaults=config)
    try:
        jobtype_list = parser.get_str("jobtype")
    except NoOptionError:
        raise KeyError("Necessary option 'jobtype' missing from config file")
    for jobtype in jobtype_list.split(","):
        if jobtype == "scan":
            from energyscan.scan import scan_main as jobtype_main
        else:
            raise WrongJobtypeError("Wrong jobtype chosen. Supported ones are: scan")
        jobtype_main(parser)

if __name__ == "__main__":
    if len(sys.argv)==1:
        _print_example()
    else:
        for infile in sys.argv[1:]:
            _main(infile)
