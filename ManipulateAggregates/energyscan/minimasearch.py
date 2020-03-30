"""Search for local energy minimum geometries.

This subsubmodule is part of ManipulateAggregates.energyscan. It implements the
second step of the 3-step procedure that creates low energy aggregate geometries.

Parallelization is supported for this subsubmodule.
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
import os, re, sys, io
from multiprocessing import Pool, Event

import logging

logger = logging.getLogger(__name__)
try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from FireDeamon import (
        RegularNeighbourListPy,
        IrregularNeighbourListPy,
        LocalMinimaPy,
    )
except ImportError:
    logger.warning(
        "Could not import RegularNeighbourListPy, IrregularNeighbourListPy or LocalMinimaPy from FireDeamon"
    )
try:
    from ..energyscan.ansilliary import general_grid, get_old_dxfiles, init_hashing
except ImportError:
    logger.warning(
        "Could not import general_grid, get_old_dxfiles or init_hashing from ..energyscan.ansilliary"
    )
try:
    from ..collection.read import read_dx
except ImportError:
    logger.warning("Could not import ..collection.read.read_dx")
try:
    from ..collection.write import CommentError
except ImportError:
    logger.warning("Could not import .collection.write.CommentError")
try:
    from ..collection import hashIO
except ImportError:
    logger.warning("Could not import ..collection.hashIO")
try:
    from ..collection.p2p3IO import open, close, writeto, hashstring
except ImportError:
    logger.warning("Could not import p2p3IO")

# default process name
PROCNAME = "EScan.MS"
try:
    from FireDeamon import set_procname
except ImportError:
    set_procname = lambda s: None
    logger.warning("Could not import FireDeamon.set_procname")

global data_ms


def _minimasearch_parallel_init(c_neighbour_list, terminating, PROCNAME):
    """Allow easy data sharing between processes without pickling.

    Args:
        c_neighbour_list: (internal data structure of libFireDeamon) a data
            structure describing which point on a grid is whose neighbour
        terminating: (return value of multiprocessing.Event()) specifies 
            whether or not a parallel computation has been terminated
            prematurely or not
        PROCNAME: (string) name of the process
    """
    global data_ms
    data_ms = (c_neighbour_list, terminating, PROCNAME)


def _minimasearch_process(args):
    """Each worker process executes this function.

    Args:
        args: (list) arguments to be passed to the worker processes (via
            pickling)
    """
    global data_ms

    c_neighbour_list, terminating, PROCNAME = data_ms

    set_procname(PROCNAME + ".%d" % (os.getpid()))

    try:
        if not terminating.is_set():
            (
                single_file,
                degeneration,
                nr_neighbours,
                progress,
                upper_cutoff,
                lower_cutoff,
                depths_sort,
                gzipped,
            ) = args
            try:
                temp = read_dx(
                    single_file,
                    grid=False,
                    data=True,
                    silent=True,
                    comments=True,
                    gzipped=gzipped,
                )
            except ValueError as e:
                print(
                    "Error when reading in dx-file %s, skipping. Error was:"
                    % (single_file),
                    e,
                    file=sys.stderr,
                )
                return None

            a1, a2, a3 = list(
                map(float, re.split(r",|\(|\)", temp["comments"][0])[1:4])
            )
            tempvalues = temp["data"]

            depths = []
            minima = LocalMinimaPy(
                c_neighbour_list,
                tempvalues,
                degeneration,
                nr_neighbours,
                prog_report=(progress == 1),
                upper_cutoff=upper_cutoff,
                lower_cutoff=lower_cutoff,
                sort_it=depths_sort,
                depths=depths,
            )
    except KeyboardInterrupt:
        print(
            "Terminating worker process " + str(os.getpid()) + " prematurely.",
            file=sys.stderr,
        )

    if depths_sort == 0:
        depths = [0.0] * len(minima)
    return minima, depths, [tempvalues[m] for m in minima], (a1, a2, a3)


def minimasearch_main(parser):
    """Main control function for the minima search procedure.

    Args:
        parser: (of class ManipulateAggregates.collection.read.SectionlessConfigParser)
            contains information about the config file. Defines the methods
            "get_str", "get_int", "get_float" and "get_boolean" to get the
            appropriate data type.
    """
    set_procname(PROCNAME)
    gets = parser.get_str
    geti = parser.get_int
    getf = parser.get_float
    getb = parser.get_boolean
    do_calculate = not (getb("config_check"))
    # do some error checking
    # value for progress reports
    if geti("progress") not in [0, 1, 2]:
        raise ValueError(
            'Wrong value for parameter "progress" given. Must be 0,1 or 2.'
        )
    else:
        progress = geti("progress")
    # check whether partitioning over nodes was switched on
    if not gets("partition") == "1/1":
        raise ValueError("Parallelization unequal 1/1 not suported for minima search.")
    # boolean values
    # NONE YET PRESENT FOR THIS JOBTYPE
    # string values
    if not gets("neighbour_check_type") in [
        "eukledian",
        "manhattan_single",
        "manhattan_multiple",
    ]:
        raise ValueError(
            "Option neighbour_check must be 'eukledian', 'manhattan_single' or 'manhattan_multiple'."
        )
    # float values (or lists of floats)
    cutoff_scale = getf("cutoff_scale")

    init_hashing(geti("hashdepth"), geti("hashwidth"), gets("hashalg"))

    for check_option in ["degeneration", "maxval", "depths_sort"]:
        getf(check_option)
    # check whether some options conflict
    # NO CONFLICTS KNOWN YET

    # spatial grid: check gridtype and set-up grid
    # read in parameters that are required in any case for the appropriate gridtypes
    if gets("sp_gridtype") in ("full", "half"):
        np_counts = numpy.array(list(map(int, gets("countsxyz").split(","))), dtype=int)
        np_del = numpy.array(list(map(float, gets("distxyz").split(","))), dtype=float)
        np_org = numpy.array([0, 0, 0], dtype=float)
    # treat auto-adjustment
    if not gets("sp_autoadjust") in ("", "none"):
        gfdict = {"TYPE": "%s" % (gets("sp_gridtype"))}
        # For each type of grid, define which parameters are allowed to be auto-adjusted
        # and how many can be adjusted at the same time.
        # Also, generate string representations of what you expect to find in the
        # gridfile in order to find out which parameters deviate.
        if gets("sp_gridtype") in ("full", "half"):
            gfdict["COUNTS"] = "%d,%d,%d" % tuple(np_counts)
            gfdict["DIST"] = "%.3f,%.3f,%.3f" % tuple(np_del)
            gfdict["ORG"] = "%.3f,%.3f,%.3f" % tuple(np_org)
            gfallowed = ("COUNTS", "DIST")
            gfmaxadjust = 1
        else:
            raise ValueError(
                "WARNING: grid auto-adjustment not supported for current gridtype %s"
                % (gets("sp_gridtype"))
            )
        if gets("sp_gridtype") == "half":
            gfdict["TYPE"] += gets("halfspace")
        # Check whether file containing information about auto-adjustment exists and
        # whether too many deviations are found.
        try:
            with open(gets("sp_gridsave"), "r") as gf:
                deviations = 0
                devstring = ""
                for line in gf:
                    l = line.rstrip().split()
                    if len(l) != 2:
                        raise ValueError(
                            "Gridtype file %s must not contain more or less than 2 columns per line."
                        )
                    if l[0] in gfallowed:
                        if gfdict[l[0]] != l[1]:
                            deviations += 1
                            devstring += l[0] + " "
                            gfdict[l[0]] = l[1]
                    else:
                        if gfdict[l[0]] != l[1]:
                            raise ValueError(
                                "Mandatory parameter %s is not equal for grid defined in config file and grid defined in grid file %s."
                                % (l[0], gets("sp_gridsave"))
                            )
                if deviations > gfmaxadjust:
                    raise ValueError(
                        "The grid used for the scan (see file %s) deviates from the one defined "
                        % (gets("sp_gridsave"))
                        + "in the config file by more than %d allowed parameter(s), namely %s. It should be:"
                        % (gfmaxadjust, devstring),
                        gfdict,
                    )
                elif deviations == 0:
                    print(
                        "...successfully read in grid used for the scan (was not auto-adjusted)..."
                    )
                else:
                    print(
                        "...successfully read in grid used for the scan (auto-adjusted in %d parameters: %s)..."
                        % (deviations, devstring)
                    )
        except (IOError, OSError) as e:
            raise IOError(
                "Auto-adjustment of spatial grid requested but file containing grid information %s could not be opened."
                % (gets("sp_gridsave")),
                e,
            )
        # For each gridtype, use the deviating values.
        if gets("sp_gridtype") in ("full", "half"):
            for ds in devstring.split():
                if ds == "COUNTS":
                    np_counts = numpy.array(
                        list(map(int, gfdict[ds].split(","))), dtype=int
                    )
                elif ds == "DIST":
                    np_del = numpy.array(
                        list(map(float, gfdict[ds].split(","))), dtype=float
                    )
    if gets("sp_gridtype") == "full":
        gets("suffix")
    elif gets("sp_gridtype") == "half":
        np_counts_pos = numpy.array([c for c in np_counts])
        np_counts_neg = numpy.array([c for c in np_counts])
        halfspace_vec = list(map(int, gets("halfspace").split(",")))
        for i in (0, 1, 2):
            if halfspace_vec[i] < 0:
                np_counts_pos[i] = abs(halfspace_vec[i])
            if halfspace_vec[i] > 0:
                np_counts_neg[i] = abs(halfspace_vec[i])
        gets("suffix")
    else:
        raise ValueError("Wrong value for config value sp_gridtype.")

    if gets("distance_cutoff") == "auto":
        if (
            gets("neighbour_check_type") == "eukledian"
            or gets("neighbour_check_type") == "manhattan_single"
        ):
            raise ValueError(
                "Value 'auto' for option 'distance_cutoff' only supported for 'manhattan_multiple' and grids 'full' or 'half'."
            )
        elif gets("neighbour_check_type") == "manhattan_multiple":
            if gets("sp_gridtype") in ("full", "half"):
                distance_cutoff = list(cutoff_scale * np_del)
            else:
                raise ValueError(
                    "Value 'auto' for option 'distance_cutoff' only supported for grids 'full' or 'half'."
                )
        else:
            raise Exception("Wrong value for option 'neighbour_check_type'.")
    else:
        if (
            gets("neighbour_check_type") == "eukledian"
            or gets("neighbour_check_type") == "manhattan_single"
        ):
            distance_cutoff = getf("distance_cutoff")
        elif gets("neighbour_check_type") == "manhattan_multiple":
            try:
                distance_cutoff = list(
                    map(
                        lambda s: cutoff_scale * float(s),
                        gets("distance_cutoff").split(","),
                    )
                )
            except ValueError:
                raise TypeError(
                    "Each element of option distance_cutoff must be of type float."
                )
            if len(distance_cutoff) != 3:
                raise ValueError(
                    "Option 'distance_cutoff' must have three entries for 'manhattan_multiple'."
                )
        else:
            raise Exception("Wrong value for option 'neighbour_check_type'.")

    # get number of neighbours to search
    nr_shells = None
    if gets("nr_neighbours") == "auto":
        if gets("neighbour_check_type") == "manhattan_multiple" and gets(
            "sp_gridtype"
        ) in ("full", "half"):
            nr_neighbours = [
                2 * int(1.0 * distance_cutoff[i] / np_del[i]) + 1 for i in range(3)
            ]
            tmp = nr_neighbours[0]
            if all((i == tmp for i in nr_neighbours)):
                nr_shells = int((tmp - 1) // 2)
            nr_neighbours = nr_neighbours[0] * nr_neighbours[1] * nr_neighbours[2] - 1
        else:
            raise ValueError(
                "Value 'auto' for 'nr_neighbours' only supported for 'manhattan_multiple' and the sp_gridtypes 'full' and 'half'."
            )
    else:
        nr_neighbours = geti("nr_neighbours")
    if gets("max_nr_neighbours") == "auto":
        max_nr_neighbours = nr_neighbours
    else:
        max_nr_neighbours = geti("max_nr_neighbours")

    if gets("volumetric_data").startswith("from_scan,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 2:
            # angular grid: check gridtype and set-up grid
            if gets("ang_gridtype") == "full":
                # these are the counts and distances for rotation
                countsposmain = numpy.array(
                    list(map(int, gets("countspos").split(","))), dtype=int
                )
                countsnegmain = numpy.array(
                    list(map(int, gets("countsneg").split(","))), dtype=int
                )
                totcounts = countsposmain + countsnegmain + 1
                nr_dx_files = 1
                for c in totcounts:
                    nr_dx_files *= c
            else:
                raise ValueError(
                    "Option 'volumetric_data' of 'from_scan' only supported for 'ang_gridtype'=='full'"
                )
            filenames = []
            reuse_ids = {}
            for f in range(1, nr_dx_files + 1):
                reuse_ids[f] = True
            # reuse_ids = {f:True for f in range(1,nr_dx_files+1)}
            if len(gets("scan_restartdirs")) > 0:
                print("...checking which dx-files are present in old directories...")
                olddxfiles = get_old_dxfiles(
                    gets("scan_restartdirs").split(","), gets("suffix")
                )
                for c in olddxfiles:
                    reuse_ids[c] = False
                # reuse_ids.update({c:False for c in olddxfiles})
                filenames += [olddxfiles[c] for c in olddxfiles]
            discard, directory = config_data
            if os.path.isdir(directory):
                print(
                    "...checking which dx-files are present in directory created by a previous scan..."
                )
                filenames += [
                    hashIO.hashpath(directory + os.sep + str(f) + "_" + gets("suffix"))
                    for f in reuse_ids
                    if reuse_ids[f]
                ]
            else:
                print(
                    "WARNING: directory %s that should contain old dx-files does not exist."
                    % (directory),
                    file=sys.stderr,
                )
            print("...determining which files are missing...")
            dx_files = sorted(
                [f for f in filenames if os.path.exists(f) and os.stat(f).st_size > 0],
                key=str.lower,
            )
            missing_dx_files = sorted(
                [
                    f.split(os.sep)[-1]
                    for f in filenames
                    if not os.path.exists(f) or os.stat(f).st_size <= 0
                ],
                key=str.lower,
            )
            if len(dx_files) != nr_dx_files:
                print(
                    "WARNING: some files that were expected to be generated by the scan in directory '%s' are missing"
                    % (directory),
                    file=sys.stderr,
                )
                if len(gets("scan_restartdirs")) > 0:
                    print(
                        "         and could not be supplied from previous runs in the directories: %s"
                        % (gets("scan_restartdirs")),
                        file=sys.stderr,
                    )
                print(
                    "Missing files: "
                    + ", ".join(
                        sorted(
                            missing_dx_files,
                            key=lambda e:
                            # e.split(os.sep)[-1] is the name of the file (NAME)
                            # NAME.split("_"+gets("suffix"))[0] is the numer of the dx-file
                            int(e.split(os.sep)[-1].split("_" + gets("suffix"))[0]),
                        )
                    ),
                    file=sys.stderr,
                )
        else:
            raise ValueError(
                'Wrong format for parameter "volumetric_data" starting with "from_scan" given. Must be "from_scan,DIR".'
            )
    elif gets("volumetric_data").startswith("dir_regex,"):
        config_data = gets("volumetric_data").split(",")
        if len(config_data) == 3:
            discard, directory, regex = config_data
            if os.path.isdir(directory):
                dx_files = [
                    f
                    for f in hashIO.listfiles(
                        directory, regex, nullsize=False, nulldepth=False
                    )
                ]
            else:
                raise ValueError(
                    "Given directory " + directory + " is not a directory."
                )
        else:
            raise ValueError(
                'Wrong format for parameter "volumetric_data" starting with "dir_regex" given. Must be "dir_regex,DIR,REGEX".'
            )
    else:
        raise ValueError(
            'Wrong value for parameter "volumetric_data" given. Must be "from_scan,DIR" or "dir_regex,DIR,REGEX".'
        )

    use_regular = (
        (nr_shells is not None)
        and (gets("sp_gridtype") in ("full", "half"))
        and (gets("neighbour_check_type") == "manhattan_multiple")
        and (gets("nr_neighbours") == "auto")
    )

    if use_regular:
        if gets("max_nr_neighbours") != "auto":
            print(
                "WARNING: the value for 'max_nr_neighbours' is not used for regular grids.",
                file=sys.stderr,
            )
        print(
            "...using fast neighbour-search algorithm for regular grid with %d neighbour shells..."
            % (nr_shells)
        )
    else:
        print("...using slow neighbour-search algorithm for irregular grid...")
        print("Conditions not fulfilled for fast algorithm:")
        regular_dict = {
            "sp_gridtype in ('full','half')": gets("sp_gridtype") in ("full", "half"),
            "neighbour_check_type == manhattan_multiple": gets("neighbour_check_type")
            == "manhattan_multiple",
            "nr_neighbours == auto": gets("nr_neighbours") == "auto",
        }
        tmpstring = "int(distance_cutoff[i]/distxyz[i]) not the same for i in {0,1,2}"
        if (
            gets("sp_gridtype") in ("full", "half")
            and gets("neighbour_check_type") == "manhattan_multiple"
            and gets("nr_neighbours") == "auto"
        ):
            regular_dict[tmpstring] = nr_shells is not None
        else:
            regular_dict[tmpstring] = True
        for reason in regular_dict:
            if not regular_dict[reason]:
                print("     %s" % (reason))
    geti("pool_chunksize")

    if not do_calculate:
        if len(dx_files) == 0:
            print(
                "WARNING: some of the dx-files you requested are currently not present, which might.",
                file=sys.stderr,
            )
            print(
                "         not be a problem since this is only a config check.",
                file=sys.stderr,
            )
        return

    if len(dx_files) == 0:
        raise RuntimeError(
            "Could not find any non-empty dx-files matching the given criteria."
        )

    # sorting the dx-files by name to always get a sorted minima file
    dx_files = sorted(
        dx_files,
        key=lambda e:
        # e.split(os.sep)[-1] is the name of the file (NAME)
        # NAME.split("_"+gets("suffix"))[0] is the numer of the dx-file
        int(e.split(os.sep)[-1].split("_" + gets("suffix"))[0]),
    )
    print("...sorted list of dx-files...")

    if gets("sp_gridtype") == "full":
        np_grid = general_grid(np_org, np_counts, np_counts, np_del)
    elif gets("sp_gridtype") == "half":
        np_grid = general_grid(np_org, np_counts_pos, np_counts_neg, np_del)
    else:
        raise ValueError("Wrong value for config value sp_gridtype.")

    pos_from_index = lambda index: np_grid[index]

    if use_regular:
        if gets("sp_gridtype") == "full":
            FD_nr_points = list(map(lambda c: 2 * c + 1, np_counts))
        elif gets("sp_gridtype") == "half":
            FD_nr_points = list(np_counts_pos + np_counts_neg + 1)
        else:
            raise ValueError("Wrong value for config value sp_gridtype.")
        c_neighbour_list = RegularNeighbourListPy(
            FD_nr_points, int(nr_shells), prog_report=False, exclude_border=True
        )
    else:
        c_neighbour_list = IrregularNeighbourListPy(
            np_grid,
            nr_neighbours,
            distance_cutoff,
            max_nr_neighbours=max_nr_neighbours,
            prog_report=(progress == 1),
            cutoff_type=gets("neighbour_check_type"),
            sort_it=False,
        )
    print("...generated neighbour list...")

    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
        nr_threads = 1
    except ValueError:
        nr_threads = 1

    # how to properly handle keyboard interrupts when multi processing has been taken from:
    # http://stackoverflow.com/questions/14579474/multiprocessing-pool-spawning-new-childern-after-terminate-on-linux-python2-7
    terminating = Event()

    # global data_ms                            #DEBUG
    # data_ms = (c_neighbour_list, terminating) #DEBUG

    args = [
        [
            single_file,
            getf("degeneration"),
            nr_neighbours,
            progress,
            getf("maxval"),
            None,
            geti("depths_sort"),
            getb("gzipped"),
        ]
        for single_file in dx_files
    ]

    if not gets("minima_file_save").endswith(".gz"):
        minima_file = io.open(gets("minima_file_save"), "wb")
    else:
        try:
            from subprocess import Popen, PIPE

            gzipprocess = Popen(
                ["gzip", "-6", "-c", "-"],
                stdin=PIPE,
                stdout=io.open(gets("minima_file_save"), "wb"),
                bufsize=4096,
            )
            minima_file = gzipprocess.stdin
        except ImportError:
            print(
                "WARNING: cannot import gzip module, will treat %s as a non-gzipped one."
                % (gets("minima_file_save")[0:-3]),
                file=sys.stderr,
            )
            minima_file = io.open(gets("minima_file_save")[0:-3], "wb")
        except OSError:
            print(
                "WARNING: cannot import gzip module, will treat %s as a non-gzipped one."
                % (filename),
                file=sys.stderr,
            )
            minima_file = io.open(gets("minima_file_save")[0:-3], "wb")

    minima_file.write(
        hashstring("#%s FF: %s\n" % (gets("minima_file_save"), gets("forcefield")))
    )

    dx_file_count = 0
    dx_file_max = len(args)

    try:
        chunksize = geti("pool_chunksize")
        while dx_file_count < dx_file_max:
            pool = Pool(
                nr_threads,
                initializer=_minimasearch_parallel_init,
                initargs=(c_neighbour_list, terminating, PROCNAME),
            )  # NODEBUG
            chunkstart = dx_file_count
            if chunkstart + chunksize > dx_file_max:
                chunkend = dx_file_max
            else:
                chunkend = chunkstart + chunksize
            # loop over all dx-files via worker processes
            for temp in pool.imap(
                _minimasearch_process, args[chunkstart:chunkend]
            ):  # NODEBUG
                # for arg in args:                                            #DEBUG
                dx_file_count += 1
                # temp = _minimasearch_process(arg)                           #DEBUG
                if temp is None:
                    if progress > 0:
                        print(
                            "Skipping dx-file %d of %d: read error"
                            % (dx_file_count, dx_file_max)
                        )
                    continue
                minima, depths, min_energies, (a1, a2, a3) = temp
                tmplen = list(map(len, (minima, depths, min_energies)))
                if min(tmplen) <= 0:
                    if progress > 0:
                        print(
                            "Skipping dx-file %d of %d: no minima found"
                            % (dx_file_count, dx_file_max)
                        )
                    continue
                if not (min(tmplen) == max(tmplen)):
                    if progress > 0:
                        print(
                            "Error while processing dx-file %d of %d"
                            % (dx_file_count, dx_file_max),
                            file=sys.stderr,
                        )
                    raise RuntimeError(
                        "Error while processing dx-file %d of %d: lists do not have equal lengths"
                        % (dx_file_count, dx_file_max)
                    )
                if progress > 0:
                    print(
                        "Processing dx-file %d of %d: %d minima"
                        % (dx_file_count, dx_file_max, len(minima))
                    )

                for minimum, depth, min_energy in zip(minima, depths, min_energies):
                    min_pos = pos_from_index(minimum)
                    minima_file.write(
                        hashstring(
                            "%10d     %15.8f %15.8f %15.8f     %15.8f %15.8f %15.8f     %15.8E   %E \n"
                            % (
                                minimum,
                                min_pos[0],
                                min_pos[1],
                                min_pos[2],
                                a1,
                                a2,
                                a3,
                                min_energy,
                                depth,
                            )
                        )
                    )
                del (
                    minimum,
                    depth,
                    min_energy,
                    minima,
                    depths,
                    min_energies,
                    a1,
                    a2,
                    a3,
                    temp,
                    min_pos,
                )
            pool.close()  # NODEBUG
            pool.join()  # NODEBUG
    except KeyboardInterrupt as e:
        print("Caught keyboard interrupt.", file=sys.stderr)
        pool.terminate()  # NODEBUG
        pool.join()  # NODEBUG
        print("Terminating main routine prematurely.", file=sys.stderr)
        minima_file.close()
        raise e
    minima_file.close()
    if gets("minima_file_save").endswith(".gz"):
        gzipprocess.wait()
