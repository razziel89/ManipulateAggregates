from ConfigParser import NoOptionError
import os, re, sys

import numpy as np

from FireDeamon import NeighbourListPy, LocalMinimaPy
from energyscan.scan import general_grid
from collection.read import read_dx

class DXFilesNotFoundError(Exception):
    pass

def minimasearch_main(parser):
    gets   = parser.get_str
    geti   = parser.get_int
    getf   = parser.get_float
    getb   = parser.get_boolean
    do_calculate = not(getb("config_check"))
    #do some error checking
    #value for progress reports
    if geti("progress") not in [0,1,2]:
        raise ValueError('Wrong value for parameter "progress" given. Must be 0,1 or 2.')
    else:
        progress = geti("progress")
    #boolean values
    #NONE YET PRESENT FOR THIS JOBTYPE
    #string values
    if not gets("neighbour_check_type") in ['eukledian','manhattan_single','manhattan_multiple']:
        raise TypeError("Option neighbour_check must be 'eukledian', 'manhattan_single' or 'manhattan_multiple'.")
    #float values (or lists of floats)
    try:
        cutoff_scale = getf("cutoff_scale")
    except ValueError:
        raise TypeError("Option cutoff_scale must be of type float.")
    if gets("neighbour_check_type") == 'eukledian' or gets("neighbour_check_type") == 'manhattan_single':
        try:
            distance_cutoff = getf("distance_cutoff")
        except ValueError:
            raise TypeError("Option distance_cutoff must be of type float.")
    elif gets("neighbour_check_type") == 'manhattan_multiple':
        try:
            distance_cutoff = list(map(lambda s: cutoff_scale*float(s),gets("distance_cutoff").split(",")))
        except ValueError:
            raise TypeError("Each element of option distance_cutoff must be of type float.")
        if len(distance_cutoff) != 3:
            raise ValueError("Option 'distance_cutoff' must have three entries for 'manhattan_multiple'.")
    else:
        raise Exception("UNHANDLED INTERNAL ERROR")

    for check_option in ["degeneration","maxval","depths_sort"]:
        try:
            getf(check_option)
        except ValueError:
            raise TypeError("Option "+check_option+"  must be of type float.")
    #check whether some options conflict
    #NO CONFLICTS KNOWN YET
    #populate all variables with the given values
    try:
        #spatial grid: check gridtype and set-up grid
        option = "sp_gridtype"
        if gets("sp_gridtype") == "full":
            #these are only the counts in one direction
            option="countsxyz"
            np_counts = np.array(map(int,gets(option).split(",")))
            #example: 0.35,0.5,0.5
            option="distxyz"
            np_del    = np.array(map(float,gets(option).split(",")))
            np_org    = np.array([0,0,0])
            option    = "suffix"
            if do_calculate:
                np_grid   = general_grid(np_org,np_counts,np_counts,np_del)
                dx_dict = {"filename": gets(option), "counts": list(2*np_counts+1), "org": list(np_grid[0]),
                           "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}
            else:
                gets(option)
        else:
            raise ValueError("Wrong value for config value sp_gridtype.")
    except NoOptionError:
        raise KeyError("Necessary option missing from config file: "+option)
    #get number of neighbours to search
    if gets("nr_neighbours") == "auto":
        if gets("neighbour_check_type") == "manhattan_multiple" and gets("sp_gridtype") == "full":
            import operator
            nr_neighbours = reduce(operator.mul,(2*int(1.0*distance_cutoff[i]/np_del[i])+1 for i in xrange(3))) - 1
        else:
            raise ValueError("Value 'auto' for 'nr_neighbours' only supported for 'manhattan_multiple' and the sp_gridtype 'full'.")
    else:
        try:
            nr_neighbours = geti("nr_neighbours")
        except ValueError:
            raise TypeError("Option nr_neighbours must be of type int.")
    if gets("max_nr_neighbours") == "auto":
        max_nr_neighbours = nr_neighbours
    else:
        try:
            max_nr_neighbours = geti("max_nr_neighbours")
        except ValueError:
            raise TypeError("Option max_nr_neighbours must be of type int or 'auto'.")

    if gets("volumetric_data").startswith("from_scan,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 2:
            #angular grid: check gridtype and set-up grid
            if gets("ang_gridtype") == "full":
                #these are the counts and distances for rotation
                countsposmain = np.array(map(int,gets("countspos").split(",")))
                countsnegmain = np.array(map(int,gets("countsneg").split(",")))
                nr_dx_files   = reduce(operator.mul,countsposmain+countsposmain+1)
            else:
                raise ValueError("Option 'volumetric_data' of 'from_scan' only supported for 'ang_gridtype'=='full'")
            discard,directory = config_data
            if os.path.isdir(directory):
                filenames = [directory+"/"+str(f)+"_"+gets("suffix") for f in xrange(1,nr_dx_files+1)]
                dx_files = [f for f in filenames if os.path.exists(f) and os.stat(f).st_size > 0]
                if len(dx_files) != nr_dx_files:
                    print >>sys.stderr,"WARNING: some files that were expected to be generated by the scan are missing."
                    missing_files=""
                    for f in filenames:
                        if f not in dx_files:
                            missing_files += f+", "
                    print >>sys.stderr,"Missing files: "+missing_files[:-2]
        else:
            raise ValueError('Wrong format for parameter "volumetric_data" starting with "from_scan" given. Must be "from_scan,DIR".')
    elif gets("volumetric_data").startswith("dir_regex,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 3:
            discard,directory,regex = config_data
            if os.path.isdir(directory):
                dx_files = [f for f in os.listdir() if re.search(regex, f) and os.stat(f).st_size > 0]
            else:
                raise ValueError('Given directory '+directory+' is not a directory.')
        else:
            raise ValueError('Wrong format for parameter "volumetric_data" starting with "dir_regex" given. Must be "dir_regex,DIR,REGEX".')
    else:
        raise ValueError('Wrong value for parameter "volumetric_data" given. Must be "from_scan,DIR" or "dir_regex,DIR,REGEX".')

    if not do_calculate:
        if len(dx_files) == 0:
            print >>sys.stderr,"WARNING: some of the dx-files you requested are currently not present, which might."
            print >>sys.stderr,"         not be a problem since this is only a config check."
        return

    if len(dx_files) == 0:
        raise DXFilesNotFoundError("Could not find any non-empty dx-files matching the given criteria.")

    c_neighbour_list = NeighbourListPy(np_grid, nr_neighbours, distance_cutoff, max_nr_neighbours=max_nr_neighbours,
                                       prog_report=(progress==1), cutoff_type=gets("neighbour_check_type"),
                                       sort_it=False)
    print "...generated neighbour list..."

    #loop over all dx-files
    for single_file in dx_files:
        values = read_dx(single_file,grid=False,data=True,silent=True)['data']
        return_depths = []
        LocalMinimaPy(c_neighbour_list, values, getf("degeneration"), nr_neighbours,
                      prog_report=(progress==1), upper_cutoff=getf("maxval"), lower_cutoff=None,
                      sort_it=geti("depths_sort"), depths=return_depths)
        del values
