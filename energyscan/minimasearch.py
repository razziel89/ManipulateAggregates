import os, re, sys
from multiprocessing import Pool, Event

import numpy as np

from FireDeamon import RegularNeighbourListPy, IrregularNeighbourListPy, LocalMinimaPy
from energyscan.scan import general_grid,get_old_dxfiles,_init_hashing
from collection.read import read_dx
from collection.write import CommentError
from collection import hashIO

global data_ms

#this allows for easy data sharing between processes without pickling
def minimasearch_parallel_init(c_neighbour_list, terminating):
    global data_ms
    data_ms = (c_neighbour_list, terminating)

class DXFilesNotFoundError(Exception):
    pass

def _minimasearch_process(args):
    global data_ms

    c_neighbour_list, terminating = data_ms

    try:
        if not terminating.is_set():
            single_file, degeneration, nr_neighbours, progress, upper_cutoff, lower_cutoff, depths_sort, gzipped = args
            try:
                temp = read_dx(single_file,grid=False,data=True,silent=True,comments=True,gzipped=gzipped)
            except ValueError as e:
                print >>sys.stderr,"Error when reading in dx-file %s, skipping. Error was:"%(single_file),e
                return None

            a1,a2,a3   = list(map(float,re.split(r',|\(|\)',temp["comments"][0])[1:4]))
            tempvalues = temp['data']

            depths = []
            minima = LocalMinimaPy(c_neighbour_list, tempvalues, degeneration, nr_neighbours,
                                   prog_report=(progress==1), upper_cutoff=upper_cutoff, lower_cutoff=lower_cutoff,
                                   sort_it=depths_sort, depths=depths)
    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

    return minima,depths,[tempvalues[m] for m in minima],(a1,a2,a3)

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
        raise ValueError("Option neighbour_check must be 'eukledian', 'manhattan_single' or 'manhattan_multiple'.")
    #float values (or lists of floats)
    cutoff_scale = getf("cutoff_scale")
    if gets("neighbour_check_type") == 'eukledian' or gets("neighbour_check_type") == 'manhattan_single':
        distance_cutoff = getf("distance_cutoff")
    elif gets("neighbour_check_type") == 'manhattan_multiple':
        try:
            distance_cutoff = list(map(lambda s: cutoff_scale*float(s),gets("distance_cutoff").split(",")))
        except ValueError:
            raise TypeError("Each element of option distance_cutoff must be of type float.")
        if len(distance_cutoff) != 3:
            raise ValueError("Option 'distance_cutoff' must have three entries for 'manhattan_multiple'.")
    else:
        raise Exception("UNHANDLED INTERNAL ERROR")

    _init_hashing(geti("hashdepth"),geti("hashwidth"),gets("hashalg"))


    for check_option in ["degeneration","maxval","depths_sort"]:
        getf(check_option)
    #check whether some options conflict
    #NO CONFLICTS KNOWN YET
    #populate all variables with the given values
    #spatial grid: check gridtype and set-up grid
    option = "sp_gridtype"
    if gets("sp_gridtype") == "full":
        #these are only the counts in one direction
        option="countsxyz"
        np_counts = np.array(map(int,gets(option).split(",")),dtype=int)
        #example: 0.35,0.5,0.5
        option="distxyz"
        np_del    = np.array(map(float,gets(option).split(",")))
        np_org    = np.array([0,0,0])
        gets("suffix")
    else:
        raise ValueError("Wrong value for config value sp_gridtype.")
    import operator
    #get number of neighbours to search
    nr_shells = None
    if gets("nr_neighbours") == "auto":
        if gets("neighbour_check_type") == "manhattan_multiple" and gets("sp_gridtype") == "full":
            nr_neighbours = [2*int(1.0*distance_cutoff[i]/np_del[i])+1 for i in xrange(3)]
            tmp = nr_neighbours[0]
            if all((i==tmp for i in nr_neighbours)):
                nr_shells = (tmp-1)/2
            nr_neighbours = reduce(operator.mul,nr_neighbours) - 1
        else:
            raise ValueError("Value 'auto' for 'nr_neighbours' only supported for 'manhattan_multiple' and the sp_gridtype 'full'.")
    else:
        nr_neighbours = geti("nr_neighbours")
    if gets("max_nr_neighbours") == "auto":
        max_nr_neighbours = nr_neighbours
    else:
        max_nr_neighbours = geti("max_nr_neighbours")

    if gets("volumetric_data").startswith("from_scan,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 2:
            #angular grid: check gridtype and set-up grid
            if gets("ang_gridtype") == "full":
                #these are the counts and distances for rotation
                countsposmain = np.array(map(int,gets("countspos").split(",")))
                countsnegmain = np.array(map(int,gets("countsneg").split(",")))
                nr_dx_files   = reduce(operator.mul,countsposmain+countsnegmain+1)
            else:
                raise ValueError("Option 'volumetric_data' of 'from_scan' only supported for 'ang_gridtype'=='full'")
            filenames = []
            reuse_ids = {f:True for f in xrange(1,nr_dx_files+1)}
            if len(gets("scan_restartdirs")) > 0:
                print "...checking which dx-files are present in old directories..."
                olddxfiles = get_old_dxfiles(gets("scan_restartdirs").split(","),gets("suffix"))
                reuse_ids.update({c:False for c in olddxfiles})
                filenames += [olddxfiles[c] for c in olddxfiles]
            discard,directory = config_data
            if os.path.isdir(directory):
                print "...checking which dx-files are present in directory created by a previous scan..."
                filenames += [hashIO.hashpath(directory+os.sep+str(f)+"_"+gets("suffix")) for f in reuse_ids if reuse_ids[f]]
            else:
                print >>sys.stderr,"WARNING: directory %s that should contain old dx-files does not exist."%(directory)
            print "...determining which files are missing..."
            dx_files         = sorted([f for f in filenames if os.path.exists(f) and os.stat(f).st_size > 0], key=str.lower)
            missing_dx_files = sorted([f.split(os.sep)[-1] for f in filenames if not os.path.exists(f) or os.stat(f).st_size <= 0], key=str.lower)
            if len(dx_files) != nr_dx_files:
                print >>sys.stderr, "WARNING: some files that were expected to be generated by the scan in directory '%s' are missing"%(directory)
                if len(gets("scan_restartdirs")) > 0:
                    print >>sys.stderr, "         and could not be supplied from previous runs in the directories: %s"%(gets("scan_restartdirs"))
                missing_files=", ".join(missing_dx_files)
                print >>sys.stderr,"Missing files: "+missing_files
        else:
            raise ValueError('Wrong format for parameter "volumetric_data" starting with "from_scan" given. Must be "from_scan,DIR".')
    elif gets("volumetric_data").startswith("dir_regex,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 3:
            discard,directory,regex = config_data
            if os.path.isdir(directory):
                dx_files = [f for f in hashIO.listfiles(directory,regex,nullsize=False,nulldepth=False)]
            else:
                raise ValueError('Given directory '+directory+' is not a directory.')
        else:
            raise ValueError('Wrong format for parameter "volumetric_data" starting with "dir_regex" given. Must be "dir_regex,DIR,REGEX".')
    else:
        raise ValueError('Wrong value for parameter "volumetric_data" given. Must be "from_scan,DIR" or "dir_regex,DIR,REGEX".')

    use_regular = ( (nr_shells is not None) and (gets("sp_gridtype")=="full") and 
            (gets("neighbour_check_type")=="manhattan_multiple") and (gets("nr_neighbours")=="auto") )

    if use_regular:
        if gets("max_nr_neighbours") != "auto":
            print >>sys.stderr,"WARNING: the value for 'max_nr_neighbours' is not used for regular grids."
        print "...using fast neighbour-search algorithm for regular grid with %d neighbour shells..."%(nr_shells)
    else:
        print "...using slow neighbour-search algorithm for irregular grid..."
        print "Conditions not fulfilled for fast algorithm:"
        regular_dict = {
                "sp_gridtype == full"                        : gets("sp_gridtype")=="full",
                "neighbour_check_type == manhattan_multiple" : gets("neighbour_check_type")=="manhattan_multiple",
                "nr_neighbours == auto"                      : gets("nr_neighbours")=="auto"
                }
        tmpstring = "int(distance_cutoff[i]/distxyz[i]) not the same for i in {0,1,2}"
        if gets("sp_gridtype")=="full" and gets("neighbour_check_type")=="manhattan_multiple" and gets("nr_neighbours")=="auto":
            regular_dict[tmpstring] = nr_shells is not None
        else:
            regular_dict[tmpstring] = True
        for reason in regular_dict:
            if not regular_dict[reason]:
                print "     %s"%(reason)

    if not do_calculate:
        if len(dx_files) == 0:
            print >>sys.stderr,"WARNING: some of the dx-files you requested are currently not present, which might."
            print >>sys.stderr,"         not be a problem since this is only a config check."
        return

    if len(dx_files) == 0:
        raise DXFilesNotFoundError("Could not find any non-empty dx-files matching the given criteria.")

    np_grid        = general_grid(np_org,np_counts,np_counts,np_del)
    pos_from_index = lambda index: np_grid[index]

    if use_regular:
        c_neighbour_list = RegularNeighbourListPy(list(map(lambda c: 2*c+1, np_counts)), int(nr_shells), prog_report=False)
    else:
        c_neighbour_list = IrregularNeighbourListPy(np_grid, nr_neighbours, distance_cutoff, max_nr_neighbours=max_nr_neighbours,
                                           prog_report=(progress==1), cutoff_type=gets("neighbour_check_type"),
                                           sort_it=False)
    print "...generated neighbour list..."

    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
        nr_threads = 1
    except ValueError:
        nr_threads = 1

    #how to properly handle keyboard interrupts when multi processing has been taken from:
    #http://stackoverflow.com/questions/14579474/multiprocessing-pool-spawning-new-childern-after-terminate-on-linux-python2-7
    terminating = Event()
    
    pool = Pool(nr_threads, initializer=minimasearch_parallel_init, initargs=(c_neighbour_list, terminating))   #NODEBUG
    #global data_ms                            #DEBUG
    #data_ms = (c_neighbour_list, terminating) #DEBUG

    args = [[ single_file,getf("degeneration"),nr_neighbours,progress,getf("maxval"),None,geti("depths_sort"), getb("gzipped")]
        for single_file in dx_files]

    minima_file = open(gets("minima_file_save"),"wb")
    minima_file.write("#%s\n"%(gets("minima_file_save")))

    dx_file_count = 1
    dx_file_max   = len(dx_files)

    try:
        if progress>0:
            print "Processing dx-file %d of %d"%(dx_file_count,dx_file_max)
        #loop over all dx-files via worker processes
        for temp in pool.imap_unordered(_minimasearch_process, args):    #NODEBUG
        #for temp in map(_minimasearch_process, args):                   #DEBUG
            dx_file_count += 1
            if progress>0:
                print "Processing dx-file %d of %d"%(dx_file_count,dx_file_max)
            if temp is None:
                continue
            minima,depths,min_energies,(a1,a2,a3) = temp

            for minimum,depth,min_energy in zip(minima,depths,min_energies):
                min_pos = pos_from_index(minimum)
                minima_file.write("%10d     %15.8f %15.8f %15.8f     %15.8f %15.8f %15.8f     %15.8E   %E \n"%(
                    minimum, min_pos[0], min_pos[1], min_pos[2], a1, a2, a3, min_energy, depth
                    ))
    except KeyboardInterrupt:
        print >>sys.stderr,"Caught keyboard interrupt."
        pool.terminate()    #NODEBUG
        print >>sys.stderr,"Terminating main routine prematurely."
    finally:
        pool.join() #NODEBUG
        #pass       #DEBUG
        minima_file.close()
