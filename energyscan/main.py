#!/usr/bin/env python
import sys
import operator

from collection.read import read_config_file as rf
from collection.read import NoOptionInConfigFileError

class WrongJobtypeError(Exception):
    pass

def _print_example():
    """
    Print an example config file for an energyscan to stdout.
    """
    s="""#This is an example config file that also tries to give some explanations about what all the parameters do.

#all lines starting with # are comments and can be removed
###VALUES NEEDED BY SEVERAL JONBTYPES AND GENERAL VALUES###
#declare the jobtype. NO DEFAULT SO IT MUST BE PROVIDED. Values are: scan, minimasearch
#May be a comma-separated list of jobtypes which will then be performed in sequence
jobtype         = scan,minimasearch,similarityscreening
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
#which geometry to use for SCAN and SIMILARITYSCREENING. If geometry1 and geometry2 are declared, this value will not be used.
geometry        = aligned.xyz
#use the same geometry again, twice (general variable replacement). This causes the programme not
# to use the value for 'geometry' for the respective geometry.
geometry1       = %(geometry)s
geometry2       = %(geometry1)s
#whether or not you want to align the molecules with their center to 0,0,0 and their third/second main axis to the x/y-axis
# prior to any calculation. optional, default: True
prealign        = True
#whether or not you want to have dxfiles written in gzipped format to save disk space. This will put more load
# on the processor when reading or writing. optional, default: False
gzipped         = False
#IMPORTANT NOTICE: declare the exact same grid for a MINIMASEARCH jobtype that was used for a previous SCAN run!

###JOBTYPE SCAN###
#give I1/I2 (with I1 and I2 being integers and I2>=I1). If not equal 1/1, this job is part of a split
#job and will only perform a certain subset (i.e., partition) of the scan (e.g. 1/5 would perform the first fifth,
#2/5 the second fifth). An exception will be raised for any jobtype other than scan if this value is not 1/1.
#optional, default: 1/1
partition       = 1/1
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
#This value must be larger than any other energy value you expect (in units of the selected
#forcefield) since all filtered values will be set to this
maxval          = 1000000000
#If True, after scanning all energies, set all values that are
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
#if this is a restarted scan, declare the directories (comma separated) where all previous data
#can be found. optional, default: EMPTY
scan_restartdirs =
#to limit the number of files per directory, this programmes uses the first letters of the hash of a
#dx-file's name to put it in subdirectories. Declare the hashing algorithm to be used (get a list
#of all supported ones via python -c 'import hashlib;print hashlib.algorithms;'), string, optional ,default: md5
hashalg         = md5
#declare how many hex-digits shall be used per directory level, optional, default: 2
hashwidth       = 2
#declare how many levels of directories shall be used for the hasing process, optional, default: 2
hashdepth       = 2

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
#           This option respects the value of 'scan_restartdirs' and will also use those dx-files.
#dir_regex: take those dx-files that match the given regular expression. They will be sorted by the first integer number
#           in the name. I.e. 'dir_regex,/home/test/dir,\\\\.dx$' would match everything ending on ".dx" in "/home/test/dir".
#           Please double backslashes. The regular expression and DIR must not contain commas.
volumetric_data   = from_scan,.
#declare the file to which the data about the minima shall be saved
minima_file_save  = minima.dat

###JOBTYPE SIMILARITYSCREENING###
#from where to take the data about the minima that were found. optional, default: same as minima_file_save
minima_file_load  = %(minima_file_save)s
#how many geometries the user wants at least. Those geometries are as diverse as possible in their geometries.
# Will try to find the number closest to the given one, but you might also get fewer depending on the cutoffs
# for RMSD and energy.
nr_geometries     = 10
#the maximum RMSD-cutoff for SIMILARITYSCREENING.
rmsd_min          = 1
#starting from rmsd_max, increase the RMSD-cutoff by this value until fewer than nr_geometries were found.
# Suppose that one is called rmsd_max, return the geometries for rmsd_max minus rmsd_step.
rmsd_step         = 0.5
#only consider geometries whose energy is closer to that of the global minimum geometry for the screenig.
# A negative value switches off screening by energy. This value is given in meV (milli electron volts).
# optional, default: -100
energy_cutoff     = -100
#whether or not to use the given energy cutoff in force field units (differs with force field). optional,
# default: False
use_ff_units      = False
#declare the xyz-file to which all geometries shall be saved that passed the screening by RMSD and possibly
# energy. optional, default screened.xyz
screened_xyz      = screened.xyz
"""
    print s

def _main(input_file):
    #default config
    config = {
            "config_check"   : "False",
            "forcefield"     : "mmff94",
            "geometry1"      : "%(geometry)s",
            "geometry2"      : "%(geometry1)s",
            "sp_gridtype"    : "full",
            "cutoff"         : "100.0",
            "vdw_scale"      : "-1.0",
            "ang_gridtype"   : "full",
            "save_dx"        : "True",
            "columns"        : "3",
            "suffix"         : "out.dx",
            "save_aligned"   : "True",
            "prealign"       : "True",
            "prefix"         : "template_",
            "aligned_suffix" : ".aligned",
            "save_noopt"     : "True",
            "gzipped"        : "False",
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
            "energy_cutoff"  : "-100",
            "screened_xyz"   : "screened.xyz",
            "minima_file_save"     : "minima.dat",
            "minima_file_load"     : "%(minima_file_save)s",
            "neighbour_check_type" : "manhattan_multiple",
            "max_nr_neighbours"    : "%(nr_neighbours)s",
            "scan_restartdirs"     : "",
            "hashdepth"            : "2",
            "hashwidth"            : "2",
            "hashalg"              : "md5",
            "use_ff_units"         : "False",
            "partition"            : "1/1",
            }
    options = [o for o in config] + [
            "jobtype"        , 
            "countsxyz"      ,
            "distxyz"        ,
            "geometry"       ,
            "countspos"      ,
            "countsneg"      ,
            "dist"           ,
            "nr_geometries"  ,
            "rmsd_min"       ,
            "rmsd_step"      
            ]
    parser = rf(input_file,defaults=config)
    unknown_options = parser.check_against(options)
    if len(unknown_options)>0:
        print "WARNING: the following are unknown lines in the config file:"
        for o in unknown_options:
            print o
        print
    del unknown_options
    if parser.get_boolean("config_check"):
        print "This is a check of the config file."
    jobtype_list = parser.get_str("jobtype")
    #jobtypes have long and short names but both shall be treated the same so the following
    #is a mapping of the long and short forms to a unified form
    jobtype_dict = {
            "scan"                : "scan",
            "s"                   : "scan",
            "minimasearch"        : "minima search",
            "ms"                  : "minima search",
            "similarityscreening" : "similarity screening",
            "ss"                  : "similarity screening"
            }
    from energyscan.scan import scan_main
    from energyscan.minimasearch import minimasearch_main
    from energyscan.similarityscreening import similarityscreening_main
    functions_dict = {
            "scan"                 : scan_main,
            "minima search"        : minimasearch_main,
            "similarity screening" : similarityscreening_main
            }
    #Jobs have to be performed in a certain order to make sense. The order is given
    #in the following dictionary (starting at 0 and increasing):
    order_dict = {
            "scan"                 : 0,
            "minima search"        : 1,
            "similarity screening" : 2 
            }
    try:
        jobtype_list = sorted(  #sort the jobtypes to be performed by an increasing order parameter
                                set(    #since each job must not be performed multiple times, use a set to get the irreducible set of jobtypes
                                        [   #create a list of tuples each consisting of the desired jobtype name and its order parameter
                                            ( jobtype_dict[e.lower()] , order_dict[jobtype_dict[e.lower()]] )
                                        for e in jobtype_list.split(",")]
                                   )
                             ,key=operator.itemgetter(1))
    except KeyError as e:
        raise ValueError("Given short or long form does not match any known jobtype: %s"%(e))
    for jobtype,discard in jobtype_list:
        jobtype_main = functions_dict[jobtype]
        print "Running %s..."%(jobtype)
        jobtype_main(parser)
        print "...finished %s\n"%(jobtype)
    if parser.get_boolean("config_check"):
        print "Config file seems fine."

if __name__ == "__main__":
    if len(sys.argv)==1:
        _print_example()
    else:
        for arg in sys.argv[1:]:
            if arg == '--help':
                _print_example()
            else:
                _main(arg)
