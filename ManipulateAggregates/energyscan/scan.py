"""Systematically determine aggregate energies.

This subsubmodule is part of ManipulateAggregates.energyscan. It implements the
first step of the 3-step procedure that creates low energy aggregate geometries.

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
import sys
import os
import copy
import re
import io
import csv
import errno
from multiprocessing import Pool, Event

import logging

logger = logging.getLogger(__name__)
try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from maagbel import doubleArray, OBAggregate, OBForceField
except ImportError:
    logger.warning(
        "Could not import doubleArray, OBAggregate or OBForceField from maagbel"
    )
try:
    from maagbel import pybel
except ImportError:
    logger.warning("Could not import maagbel.pybel")
try:
    from ..aggregate import read_from_file
except ImportError:
    logger.warning("Could not import ..aggregate.read_from_file")
try:
    from .ansilliary import (
        CURSOR_UP_ONE,
        ERASE_LINE,
        init_hashing,
        double_dist,
        double_array,
    )
except ImportError:
    logger.warning(
        "Could not import CURSOR_UP_ONE, ERASE_LINE, init_hashing, double_dist or double_array from ..energyscan.ansilliary"
    )
try:
    from .ansilliary import (
        print_dx_file,
        general_grid,
        prepare_molecules,
        get_old_dxfiles,
        no_none_string,
    )
except ImportError:
    logger.warning(
        "Could not import print_dx_file, general_grid, prepare_molecules, get_old_dxfiles or no_none_string from ..energyscan.ansilliary"
    )
try:
    from ..collection.read import read_dx
except ImportError:
    logger.warning("Could not import ..collection.read.read_dx")
try:
    from ..collection import hashIO
except ImportError:
    logger.warning("Could not import ..collection.read.hashIO")
try:
    from ..collection.p2p3IO import open, writeto, close, tounicode, csvwriteopen
except ImportError:
    logger.warning("Could not import p2p3IO")

# default process name
PROCNAME = "EScan.S"
try:
    from FireDeamon import set_procname
except ImportError:
    set_procname = lambda s: None
    logger.warning("Could not import FireDeamon.set_procname")

# global data to allow easy sharing of data between processes
global data_s
data_s = None  # initialized to none
global grid

# this allows for easy data sharing between processes without pickling
def _transrot_parallel_init(obmol, transgrid, terminating, PROCNAME, mask):
    """Allow easy data sharing between processes without pickling.

    Args:
        obmol: (OpenBabel OBMol object): molecule to be treated
        transgrid: (list of C-arrays of type double) spatial grid for the scan
        terminating: (return value of multiprocessing.Event()) specifies 
            whether or not a parallel computation has been terminated
            prematurely or not
        PROCNAME: (string) name of the process
        mask: (numpy array of type bool) True or False depending on whether a
            point is masked or not
    """
    global data_s
    data_s = (obmol, transgrid, terminating, PROCNAME, mask)


def _gen_trans_en(
    obmol, obff, double_grid, maxval, cutoff, vdw_scale, report, reportstring, notmasked
):
    """Internal function, undocumented."""
    if report:
        print("   %s  %.2f%%" % (reportstring, 0.0 / len(grid)) + CURSOR_UP_ONE)
    count = 0
    transfunc = obmol.TranslatePart
    # vdwfunc    = obmol.IsGoodVDW
    vdwfunc = obmol.MinVDWDist
    setupfunc = obff.Setup
    energyfunc = obff.Energy
    for newvec, retvec in double_grid:
        if notmasked(count):
            transfunc(0, newvec)
            dist = vdwfunc(
                True, vdw_scale
            )  # True will cause maagbel to interrupt whenever a vdW clash was found
            if dist > 0.0 and dist < cutoff:
                # if vdwfunc():
                setupfunc(obmol)
                yield energyfunc(
                    False
                )  # this will cause maagbel to not evaluate gradients
            else:
                yield maxval
            transfunc(0, retvec)
        else:
            yield maxval
        count += 1
        if report and count % 1000 == 0:
            print(
                "   %s  %.2f%%" % (reportstring, 100.0 * count / len(grid))
                + CURSOR_UP_ONE
            )
    if report:
        print("   %s  %.2f%%" % (reportstring, 100.0) + CURSOR_UP_ONE)


def _trans_en(
    obmol,
    obff,
    double_grid,
    maxval,
    cutoff,
    vdw_scale,
    report=False,
    reportstring="",
    notmasked=lambda i: True,
):
    """Internal function, undocumented."""
    return list(
        _gen_trans_en(
            obmol,
            obff,
            double_grid,
            maxval,
            cutoff,
            vdw_scale,
            report,
            reportstring,
            notmasked,
        )
    )


def _transrot_en_process(args):
    """Each worker process executes this function.

    Args:
        args: (list) arguments to be passed to the worker processes (via
            pickling)
    """
    global data_s

    defaultobmol, transgrid, terminating, PROCNAME, mask = data_s

    if mask is not None:
        notmasked = lambda i: mask[i]
    else:
        notmasked = lambda i: True

    set_procname(PROCNAME + ".%d" % (os.getpid()))

    try:
        if not terminating.is_set():

            (
                a1,
                a2,
                a3,
                ffname,
                report,
                maxval,
                dx_dict,
                correct,
                savetemplate,
                templateprefix,
                anglecount,
                count,
                save_noopt,
                save_opt,
                optsteps,
                cutoff,
                vdw_scale,
                oldfile,
            ) = args

            angle_string = str(a1) + "," + str(a2) + "," + str(a3)
            angle_comment = "angles=(" + angle_string + ")"

            if oldfile is not None:
                compute = False
                try:
                    old = read_dx(
                        oldfile,
                        grid=False,
                        data=True,
                        silent=True,
                        comments=True,
                        gzipped=dx_dict["gzipped"],
                    )
                except ValueError as e:
                    print(
                        "Error when reading in old dx-file %s, recomputing. Error was:"
                        % (oldfile),
                        e,
                        file=sys.stderr,
                    )
                    compute = True
                if not compute:
                    old_a1, old_a2, old_a3 = list(
                        map(float, re.split(r",|\(|\)", old["comments"][0])[1:4])
                    )
                    if not (a1, a2, a3) == (old_a1, old_a2, old_a3):
                        print(
                            "WARNING: old dx-file %s treated %s with index %d. This is also my index but I treat %s. Recomputing."
                            % (oldfile, old["comments"][0], anglecount, angle_comment),
                            file=sys.stderr,
                        )
                        compute = True
                    else:
                        energies = old["data"].tolist()
                        del old
                if not compute:
                    if not len(transgrid) == len(energies):
                        print(
                            "WARNING: old dx-file %s contains %d entries but the spatial grid is supposed to have %d entries. Recomputing."
                            % (oldfile, len(energies), len(transgrid)),
                            file=sys.stderr,
                        )
                        compute = True
            else:
                compute = True

            if compute or savetemplate:
                obmol = OBAggregate(defaultobmol)
                obff = OBForceField.FindForceField(ffname)
                rotfunc = obmol.RotatePart
                rotfunc(0, 1, a1)
                rotfunc(0, 2, a2)
                rotfunc(0, 3, a3)

            if compute:
                energies = _trans_en(
                    obmol,
                    obff,
                    transgrid,
                    maxval * 1.2,
                    cutoff,
                    vdw_scale,
                    report=report,
                    notmasked=notmasked,
                )

                if correct or dx_dict["save_dx"]:
                    # create a copy which can then be changed and possibly saved
                    tempenergies = copy.copy(energies)

                if correct:
                    try:
                        actualmax = max((e for e in tempenergies if not e >= maxval))
                    except ValueError:
                        actualmax = maxval
                    tempenergies = [
                        actualmax if e >= maxval else e for e in tempenergies
                    ]

                if dx_dict["save_dx"]:
                    print_dx_file(
                        str(anglecount) + "_",
                        True,
                        dx_dict,
                        tempenergies,
                        angle_comment,
                    )

                if correct or dx_dict["save_dx"]:
                    del tempenergies

            if savetemplate:
                minindex = energies.index(min(energies))
                template = grid[minindex][0]

                obmol.TranslatePart(0, template)
                if obmol.IsGoodVDW(vdw_scale):
                    if save_noopt:
                        filename = (
                            templateprefix
                            + str(anglecount)
                            + "_"
                            + angle_string
                            + ".xyz"
                        )
                        pybel.Molecule(obmol).write("xyz", filename, overwrite=True)
                    if save_opt:
                        filename = (
                            templateprefix
                            + "opt_"
                            + str(anglecount)
                            + "_"
                            + angle_string
                            + ".xyz"
                        )
                        p_tempmol = pybel.Molecule(obmol)
                        p_tempmol.localopt(forcefield=ffname, steps=optsteps)
                        p_tempmol.write("xyz", filename, overwrite=True)

            # returning the molecule to its original state is not necessary since every worker process
            # creates its own instance and leaves the original one as is

    except KeyboardInterrupt:
        print(
            "Terminating worker process " + str(os.getpid()) + " prematurely.",
            file=sys.stderr,
        )

    return anglecount, (a1, a2, a3), energies, minindex


def _transrot_en(
    obmol,
    ffname,
    transgrid,
    rotgrid,
    maxval,
    dx_dict,
    correct,
    cutoff,
    vdw_scale,
    report=0,
    reportcount=1,
    reportmax=None,
    savetemplate=True,
    templateprefix="template_",
    save_noopt=True,
    save_opt=True,
    optsteps=500,
    olddxfiles={},
    partition=(1, 1),
    mask=None,
):
    """Main controller function for the actual scanning procedure."""

    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
        nr_threads = 1
    except ValueError:
        nr_threads = 1

    herereport = False
    if report == 1:
        if nr_threads == 1:
            herereport = True
        else:
            print(
                "You requested detailed progress reports but that is only supported for single threaded",
                file=sys.stderr,
            )
            print(
                "calculations (you chose "
                + str(nr_threads)
                + " threads). Will switch to semi-detailed reports.",
                file=sys.stderr,
            )

    # how to properly handle keyboard interrupts when multi processing has been taken from:
    # http://stackoverflow.com/questions/14579474/multiprocessing-pool-spawning-new-childern-after-terminate-on-linux-python2-7
    terminating = Event()

    pool = Pool(
        nr_threads,
        initializer=_transrot_parallel_init,
        initargs=(obmol, transgrid, terminating, PROCNAME, mask),
    )  # NODEBUG
    # global data_s                            #DEBUG
    # data_s = (obmol, transgrid, terminating, PROCNAME, mask) #DEBUG

    nr_angles = len(rotgrid)
    nr_points = len(transgrid)

    if len(olddxfiles) == 0:
        keyfunc = lambda c, d: d
    else:
        keyfunc = olddxfiles.get

    if partition[1] > nr_angles:
        raise ValueError(
            "Number of partitions %d cannot be greater than the number of angles %d"
            % (partition[1], nr_angles)
        )

    if not partition == (1, 1):
        print(
            "...this is a partitioned calculation (partition %d out of %d)..."
            % partition
        )

        mypartition = partition[0]
        angles_per_partition = int(nr_angles // (partition[1]))
        additional_angles = nr_angles - int(angles_per_partition * partition[1])
        partition_start_end = [None] * (partition[1])
        partition_start_end[0] = (
            0,
            angles_per_partition + (1 if additional_angles > 0 else 0),
        )
        additional_angles -= 1
        for p in range(1, partition[1]):
            partition_start_end[p] = (
                partition_start_end[p - 1][1],
                partition_start_end[p - 1][1]
                + angles_per_partition
                + (1 if additional_angles > 0 else 0),
            )
            additional_angles -= 1

        print("...partitions are:")
        for (ps, pe), c in zip(partition_start_end, range(1, partition[1] + 1)):
            print(
                "          %2d: START: %6d - END: %6d" % (c, ps + 1, pe)
                + ("   <-- that's me" if c == mypartition else "")
            )
        mypartition = (
            partition_start_end[mypartition - 1][0],
            partition_start_end[mypartition - 1][1],
        )
        reportmax = mypartition[1] - mypartition[0]
    else:
        print("...this is no partitioned calculation...")
        partition_start_end = [(0, nr_angles)]
        mypartition = (0, nr_angles)

    args = [
        [
            a1,
            a2,
            a3,
            ffname,
            herereport,
            maxval,
            dx_dict,
            correct,
            savetemplate,
            templateprefix,
            anglecount,
            count,
            save_noopt,
            save_opt,
            optsteps,
            cutoff,
            vdw_scale,
            keyfunc(anglecount, None),
        ]
        for (a1, a2, a3), anglecount, count in zip(
            rotgrid, range(reportcount, nr_angles + reportcount), range(nr_angles)
        )
    ][mypartition[0] : mypartition[1]]

    # pre-declare variables
    # the optimum energies
    opt_energies = numpy.ones((nr_points,), dtype=float) * maxval
    # the angles of the geometries corresponding to the optimum energies
    # 360 is the hard-coded default
    opt_angles = numpy.ones((nr_points, 3), dtype=float) * 360.0
    # an element is True if at the corresponding spatial point at least
    # one energy smaller than maxval has been found
    opt_present = numpy.zeros((nr_points,), dtype=bool)
    # save the optimum index for each angular arrangement
    opt_angindex = numpy.zeros((nr_angles,), dtype=int)
    # save the index of the angular arrangement that is the current optimum at the spatial point
    opt_spindex = numpy.zeros((nr_points,), dtype=int)
    # an helper array for the comparison
    np_compare = numpy.zeros((nr_points,), dtype=bool)
    anglecount = reportcount
    try:
        if mypartition[0] > 0:
            reportstring = "START --> skipping %d" % (mypartition[0])
        else:
            reportstring = "START"
        print(reportstring)
        # The structure of temp is: anglecount,(a1,a2,a3),energies,minindex
        # The function _transrot_en_process is guarantueed to return values smaller than
        # maxval only if an actual evaluation using a force field has been performed
        for temp in pool.imap_unordered(_transrot_en_process, args):  # NODEBUG
            # for arg in args:                                               #DEBUG
            # temp = _transrot_en_process(arg)                           #DEBUG
            # transform energies to numpy array
            opt_temp = numpy.array(temp[2])
            # save the optimum index of this angular arrangement for later use
            opt_angindex[temp[0] - reportcount] = temp[3]
            # get all positions where the new energies are smaller than the last optimum
            np_compare = numpy.less(opt_temp, opt_energies)
            # asign those new energies at every point where the comparison was true
            opt_energies[np_compare] = opt_temp[np_compare]
            # asign the angles at every such point
            opt_angles[np_compare] = numpy.array(temp[1])
            # find out where at least one such assignment has been performed
            opt_present[np_compare] = True
            # which angular arrangement is the current optimum at this spatial point
            # (index)
            opt_spindex[np_compare] = temp[0] - reportcount
            # result.append(temp)
            reportstring = "%d/%d == #%d" % (
                anglecount,
                reportmax,
                mypartition[0] + anglecount,
            )
            if report != 0:
                print(reportstring)
            anglecount += 1
        if nr_angles - mypartition[0] - reportmax > 0:
            reportstring = "END --> skipping %d" % (
                nr_angles - mypartition[0] - reportmax
            )
        else:
            reportstring = "END"
        print(reportstring)
        pool.close()  # NODEBUG
        pool.join()  # NODEBUG
    except KeyboardInterrupt as e:
        print("Caught keyboard interrupt.", file=sys.stderr)
        pool.terminate()  # NODEBUG
        pool.join()  # NODEBUG
        print("Terminating main routine prematurely.", file=sys.stderr)
        raise e
    return opt_energies, opt_angles, opt_spindex, opt_present, opt_angindex


def _sp_opt(
    dx,
    xyz,
    ang,
    dx_dict,
    correct,
    remove,
    maxval,
    globalopt,
    obmol,
    grid,
    transrot_result,
):
    """
    This function selects optimum geometries and energies for all
    spatial coordinates. It can also sort out such geometries that clashed
    or where the molecules were too far apart from each other.
    """
    dx_bool = no_none_string(dx)
    xyz_bool = no_none_string(xyz)
    ang_bool = no_none_string(ang)

    opt_energies, opt_angles, opt_spindex, opt_present, opt_angindex = transrot_result

    if ang_bool:
        # how to write csv files taken from https://docs.python.org/2/library/csv.html
        with csvwriteopen(ang) as csvfile:
            spamwriter = csv.writer(
                csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
            )
            if remove:
                iterable = zip(opt_angles[opt_present], grid[opt_present])
            else:
                iterable = zip(opt_angles, grid)

            for (a1, a2, a3), (x, y, z) in iterable:
                spamwriter.writerow([a1, a2, a3, x, y, z])

    if xyz_bool:
        writefile = pybel.Outputfile("xyz", xyz, overwrite=True)

        filename = xyz
        tempmol = OBAggregate(obmol)
        pybeltempmol = pybel.Molecule(tempmol)
        rotfunc = tempmol.RotatePart
        transfunc = tempmol.TranslatePart
        commentfunc = tempmol.SetTitle

        tempgrid = numpy.copy(grid)
        tempgrid[numpy.logical_not(opt_present)] = numpy.array([0.0, 0.0, 0.0])

        if remove:
            iterable = zip(
                opt_angles[opt_present],
                opt_energies[opt_present],
                tempgrid[opt_present],
            )
        else:
            iterable = zip(opt_angles, opt_energies, tempgrid)

        vec = double_array([0.0, 0.0, 0.0])
        for (a1, a2, a3), e, (x, y, z) in iterable:
            commentfunc("Energy: " + str(e))
            rotfunc(0, 1, a1)
            rotfunc(0, 2, a2)
            rotfunc(0, 3, a3)
            vec[0] = x
            vec[1] = y
            vec[2] = z
            tempmol.TranslatePart(0, vec)
            writefile.write(pybeltempmol)
            vec[0] = -x
            vec[1] = -y
            vec[2] = -z
            tempmol.TranslatePart(0, vec)
            rotfunc(0, 3, -a3)
            rotfunc(0, 2, -a2)
            rotfunc(0, 1, -a1)

        del tempmol
        writefile.close()

    if dx_bool:
        if correct:
            if len(opt_present) > 0:
                try:
                    actualmax = numpy.amax(opt_energies[opt_present])
                    values = numpy.ones(opt_energies.shape, dtype=float) * actualmax
                    values[opt_present] = opt_energies[opt_present]
                except ValueError:
                    # catch a case where no values can be assigned when
                    # opt_energies[opt_present] does not have a maximum
                    print(
                        'WARNING: cannot correct dx-files as requested by "sp_opt" and "sp_opt_dx":\n'
                        "         the maxumum value will not be the actual maximum.",
                        file=sys.stderr,
                    )
                    values = opt_energies
            else:
                print(
                    'WARNING: cannot correct dx-files as requested by "sp_opt" and "sp_opt_dx":\n'
                    "         the maxumum value will not be the actual maximum.",
                    file=sys.stderr,
                )
                values = opt_energies
        else:
            values = opt_energies
        print_dx_file(
            "", False, dx_dict, values, "Optimum energies for all spatial points."
        )

    if globalopt:
        minindex = numpy.argmin(opt_energies)
        minvalue = opt_energies[minindex]
        mina1, mina2, mina3 = opt_angles[minindex]
        minvec = double_array(grid[minindex])
        tempmol = OBAggregate(obmol)
        tempmol.RotatePart(0, 1, mina1)
        tempmol.RotatePart(0, 2, mina2)
        tempmol.RotatePart(0, 3, mina3)
        tempmol.TranslatePart(0, minvec)
        tempmol.SetTitle("Energy: " + str(minvalue))
        pybel.Molecule(tempmol).write("xyz", "globalopt.xyz", overwrite=True)
        del tempmol


def scan_main(parser):
    """Main control function for the scanning procedure.

    Args:
        parser: (of class ManipulateAggregates.collection.read.SectionlessConfigParser)
            contains information about the config file. Defines the methods
            "get_str", "get_int", "get_float" and "get_boolean" to get the
            appropriate data type.
    """
    global grid
    set_procname(PROCNAME)
    gets = parser.get_str
    geti = parser.get_int
    getf = parser.get_float
    getb = parser.get_boolean
    do_calculate = not (getb("config_check"))
    # do some error checking
    # forcefield
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError(
            'Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.'
        )
    temp_ff = OBForceField.FindType(gets("forcefield").lower())
    if temp_ff is None:
        raise ValueError(
            "Somehow there was an error loading the forcefield %s (although it should be known to OpenBabel)."
            % (gets("forcefield").lower())
        )
    del temp_ff
    # boolean values
    for check_option in [
        "save_dx",
        "save_aligned",
        "save_noopt",
        "save_opt",
        "correct",
        "sp_opt",
        "sp_correct",
        "sp_remove",
        "globalopt",
        "prealign",
        "gzipped",
    ]:
        getb(check_option)
    # remaining float values
    for check_option in ["cutoff", "vdw_scale", "maxval", "cutoff", "vdw_scale"]:
        getf(check_option)
    # remaining integer values
    for check_option in ["columns", "optsteps", "progress", "hashwidth", "hashdepth"]:
        geti(check_option)
    # check whether some options conflict
    if (
        gets("volumetric_data").startswith("from_scan,")
        and "minimasearch" in gets("jobtype").split(",")
        and not getb("save_dx")
    ):
        print(
            "WARNING: a subsequent minimasearch tries to get its dx-files from this scan but",
            file=sys.stderr,
        )
        print(
            "         you requested not to save dx-files. This is probably an error (but not so if",
            file=sys.stderr,
        )
        print(
            "         you requested those dx-files to be used from a different directory) so please check.",
            file=sys.stderr,
        )

    # initialize directory name hashing
    init_hashing(geti("hashdepth"), geti("hashwidth"), gets("hashalg"))

    # value for progress reports
    if geti("progress") not in [0, 1, 2]:
        raise ValueError(
            'Wrong value for parameter "progress" given. Must be 0,1 or 2.'
        )

    # populate all variables with the given values
    # read in the two molecules/aggregates from the given files
    mol1 = read_from_file(gets("geometry1"), ff=None)
    mol2 = read_from_file(gets("geometry2"), ff=None)

    # Compute radii of spheres that completely encompass both molecules to be able to
    # auto-adjust the gridsize
    mol1_vdw = mol1.get_vdw_radii()
    mol2_vdw = mol2.get_vdw_radii()
    mol1_coords = mol1.get_coordinates()
    mol2_coords = mol2.get_coordinates()
    mol1_center = numpy.mean(numpy.array(mol1_coords, dtype=float), axis=0)
    mol2_center = numpy.mean(numpy.array(mol2_coords, dtype=float), axis=0)
    # Will contain the radius of a sphere centered at the molecular center that
    # completely encompasses mol1
    maxdist1 = 0.0
    # Will contain the radius of a sphere centered at the molecular center that
    # completely encompasses mol2
    maxdist2 = 0.0
    distcutoff = getf("cutoff")
    vdwscale = getf("vdw_scale")
    for c1, vdw1 in zip(mol1_coords, mol1_vdw):
        npc1 = numpy.array(c1, dtype=float)
        if numpy.linalg.norm(npc1 - mol1_center) + (vdw1 * vdwscale) > maxdist1:
            maxdist1 = numpy.linalg.norm(npc1 - mol1_center) + (vdw1 * vdwscale)
    for c2, vdw2 in zip(mol2_coords, mol2_vdw):
        npc2 = numpy.array(c2, dtype=float)
        if numpy.linalg.norm(npc2 - mol2_center) + (vdw2 * vdwscale) > maxdist2:
            maxdist2 = numpy.linalg.norm(npc2 - mol2_center) + (vdw2 * vdwscale)
    # The radius of the sphere outside which no points need to be considered
    # with reduced precision and rounded up
    bigspherestr = "%.3f" % (maxdist1 + maxdist2 + distcutoff + 0.0005)
    bigsphererad = float(bigspherestr)

    # treat grid auto adjustments
    if gets("sp_gridtype") in ("full", "half"):
        # these are only the counts in one direction
        np_counts = numpy.array(list(map(int, gets("countsxyz").split(","))), dtype=int)
        np_del = numpy.array(list(map(float, gets("distxyz").split(","))), dtype=float)
        np_org = numpy.array([0, 0, 0], dtype=float)
        newdistxyz = numpy.array(
            list(
                map(
                    lambda f: float("%.3f" % (f)),
                    (bigsphererad + 0.0005) / (np_counts - 1),
                )
            ),
            dtype=float,
        )
        newcountsxyz = (
            numpy.array(list(map(int, (bigsphererad + 0.0005) / np_del)), dtype=int) + 1
        )
        distsdiffer = numpy.linalg.norm(newdistxyz - np_del) > 0.0
        countsdiffer = numpy.linalg.norm(newcountsxyz - np_counts) > 0.0
        tmpstring = "..."
        if gets("sp_autoadjust") == "distxyz":
            if distsdiffer:
                tmpstring += "changing grid size in some directions by adjusting grid spacing 'distxyz' by: "
                tmpstring += "x: %+.3f, y: %+.3f, z: %+.3f" % (
                    tuple(newdistxyz - np_del)
                )
                np_del = newdistxyz
            else:
                tmpstring += "grid dimensions need no adjustment"
        elif gets("sp_autoadjust") == "countsxyz":
            if countsdiffer:
                tmpstring += "changing grid size in some directions by adjusting number of points 'countsxyz' by: "
                tmpstring += "x: %+d, y: %+d, z: %+d" % (
                    tuple(newcountsxyz - np_counts)
                )
                np_counts = newcountsxyz
            else:
                tmpstring += "grid dimensions need no adjustment"
        elif gets("sp_autoadjust") in ("", "none"):
            if countsdiffer or distsdiffer:
                tmpstring += "won't adjust, but grid dimensions inappropriate "
                tmpstring += "in some directions by (distxyz/countsxyz): "
                tmpstring += "x: %+.3f/%+d, y: %+.3f/%+d, z: %+.3f/%+d" % (
                    tuple(
                        a
                        for b in zip(newdistxyz - np_del, newcountsxyz - np_counts)
                        for a in b
                    )
                )
            else:
                tmpstring += "grid dimensions need no adjustment"
        else:
            raise ValueError("Wrong value for config value sp_autoadjust.")
        tmpstring += "..."
        print(tmpstring)
        with open(gets("sp_gridsave"), "w") as gf:
            writeto(gf, "TYPE %s" % (gets("sp_gridtype")))
            if gets("sp_gridtype") == "half":
                writeto(gf, "%s" % (gets("halfspace")))
            writeto(gf, "\n")
            writeto(gf, "COUNTS %d,%d,%d\n" % tuple(np_counts))
            writeto(gf, "DIST %.3f,%.3f,%.3f\n" % tuple(np_del))
            writeto(gf, "ORG %.3f,%.3f,%.3f\n" % tuple(np_org))
    else:
        if not gets("sp_autoadjust") in ("", "none"):
            print(
                "WARNING: grid auto-adjustment not supported for current gridtype %s"
                % (gets("sp_gridtype")),
                file=sys.stderr,
            )

    # spatial grid: check gridtype and set-up grid
    if gets("sp_gridtype") == "full":
        if do_calculate:
            np_grid = general_grid(np_org, np_counts, np_counts, np_del)
            dx_dict = {
                "filename": gets("suffix"),
                "counts": list(2 * np_counts + 1),
                "org": list(np_grid[0]),
                "delx": [np_del[0], 0.0, 0.0],
                "dely": [0.0, np_del[1], 0.0],
                "delz": [0.0, 0.0, np_del[2]],
            }
            dx_dict["save_dx"] = getb("save_dx")
            dx_dict["gzipped"] = getb("gzipped")
        else:
            gets("suffix")
            getb("save_dx")
    elif gets("sp_gridtype") == "half":
        np_counts_pos = numpy.array([c for c in np_counts])
        np_counts_neg = numpy.array([c for c in np_counts])
        halfspace_vec = list(map(int, gets("halfspace").split(",")))
        for i in (0, 1, 2):
            if halfspace_vec[i] < 0:
                np_counts_pos[i] = abs(halfspace_vec[i])
            if halfspace_vec[i] > 0:
                np_counts_neg[i] = abs(halfspace_vec[i])
        if do_calculate:
            np_grid = general_grid(np_org, np_counts_pos, np_counts_neg, np_del)
            dx_dict = {
                "filename": gets("suffix"),
                "counts": list(np_counts_pos + np_counts_neg + 1),
                "org": list(np_grid[0]),
                "delx": [np_del[0], 0.0, 0.0],
                "dely": [0.0, np_del[1], 0.0],
                "delz": [0.0, 0.0, np_del[2]],
            }
            dx_dict["save_dx"] = getb("save_dx")
            dx_dict["gzipped"] = getb("gzipped")
        else:
            gets("suffix")
            getb("save_dx")
    else:
        raise ValueError("Wrong value for config value sp_gridtype.")
    # check whether this gives an error
    restarted = len(gets("scan_restartdirs")) > 0
    if restarted:
        olddirs = gets("scan_restartdirs").split(",")
        for d in olddirs:
            if not os.path.isdir(d):
                if do_calculate:
                    print(
                        "WARNING: directory supposed to contain dx files from previous runs %s does not exist. Skipping."
                        % (d),
                        file=sys.stderr,
                    )
                else:
                    raise ValueError(
                        "Directory supposed to contain dx files from previous runs %s does not exist."
                        % (d)
                    )
    # angular grid: check gridtype and set-up grid
    if gets("ang_gridtype") == "full":
        # these are the counts and distances for rotation
        countsposmain = numpy.array(
            list(map(int, gets("countspos").split(","))), dtype=int
        )
        countsnegmain = numpy.array(
            list(map(int, gets("countsneg").split(","))), dtype=int
        )
        distmain = numpy.array(list(map(float, gets("dist").split(","))), dtype=float)
        if do_calculate:
            np_rot = general_grid(
                numpy.array([0.0, 0.0, 0.0]), countsposmain, countsnegmain, distmain
            )
    else:
        raise ValueError("Wrong value for config value ang_gridtype.")

    partition = tuple(map(int, gets("partition").split("/")))
    if len(partition) != 2:
        raise ValueError(
            "Format for 'partition' must be I1/I2 with I1 and I2 positive integers and I1<=I2"
        )
    if partition[0] > partition[1] or partition[0] < 1 or partition[1] < 1:
        raise ValueError(
            "Format for 'partition' must be I1/I2 with I1 and I2 positive integers and I1<=I2"
        )

    # Create a mask for points whose energy never has to be evaluated
    # A "False" associated with a point means "do not evaluate its energy".
    try:
        numpy.sqrt(
            numpy.einsum(
                "ij,ij->i",
                numpy.array([[1.0, 1.0, 0.0]]),
                numpy.array([[1.0, 1.0, 0.0]]),
            )
        )
        numpy.linalg.norm(numpy.array([[1.0, 1.0, 0.0]]), axis=1)
    except AttributeError:
        normaxisone = lambda array: numpy.apply_along_axis(numpy.linalg.norm, 1, array)
    except TypeError:
        normaxisone = lambda array: numpy.sqrt(numpy.einsum("ij,ij->i", array, array))
    else:
        normaxisone = lambda array: numpy.linalg.norm(array, axis=1)
    mask = numpy.ones((len(np_grid),), dtype=bool)
    dist = numpy.zeros((len(np_grid),), dtype=float)
    origin = numpy.array([0.0, 0.0, 0.0], dtype=float)
    min1_vdw = min(mol1_vdw)
    min2_vdw = min(mol2_vdw)
    for c1, vdw1 in zip(mol1_coords, mol1_vdw):
        npc1 = numpy.array(c1, dtype=float)
        vdw = (vdw1 + min2_vdw) * vdwscale
        dist = normaxisone(np_grid - (npc1 - mol1_center))
        # dist = numpy.linalg.norm((np_grid-(npc1-mol1_center)),axis=1)
        mask[dist < vdw] = False
    inmasked = numpy.sum(mask == False)
    dist = normaxisone(np_grid - origin)
    # dist      = numpy.linalg.norm((np_grid-origin),axis=1)
    mask[dist > bigsphererad] = False
    outmasked = numpy.sum(mask == False) - inmasked
    print(
        "...computed mask, reduction in points: %.2f%%, inside: %.2f%%, outside: %.2f%%..."
        % (
            100.0 * (inmasked + outmasked) / len(mask),
            100.0 * inmasked / len(mask),
            100.0 * outmasked / len(mask),
        )
    )

    if not do_calculate:
        return

    # align the two molecules and append one to the other
    # after this, mol1 and mol2 can no longer be used
    obmol = prepare_molecules(
        mol1,
        mol2,
        gets("aligned_suffix"),
        save_aligned=getb("save_aligned"),
        align=getb("prealign"),
    )

    # convert the grid to C data types
    grid = double_dist(np_grid)

    if restarted:
        print(
            "This is a restarted run (old files are in: %s)"
            % (gets("scan_restartdirs"))
        )
        olddxfiles = get_old_dxfiles(
            gets("scan_restartdirs").split(","), gets("suffix")
        )
        print("Number of already existing dx files: %d" % (len(olddxfiles)))
    else:
        olddxfiles = {}

    # For every angle, scan the entire spatial grid and save
    # each optimum geometry if desired
    # Will also return a structure making it easy to find the optimum
    # for every spatial point
    transrot_result = _transrot_en(
        obmol,
        gets("forcefield").lower(),
        grid,
        np_rot,
        getf("maxval"),
        dx_dict,
        getb("correct"),
        getf("cutoff"),
        getf("vdw_scale"),
        report=geti("progress"),
        reportmax=len(np_rot),
        save_noopt=getb("save_noopt"),
        save_opt=getb("save_opt"),
        optsteps=geti("optsteps"),
        olddxfiles=olddxfiles,
        partition=partition,
        mask=mask,
    )

    del grid  # the grid in C data types is no longer needed since the scan has already been performed

    # Evaluate transrot_result to find the angular optimum for every
    # spatial grid point, if so desired
    if getb("sp_opt"):

        dx_dict["filename"] = gets("sp_opt_dx")
        dx_dict["save_dx"] = getb("sp_opt")

        _sp_opt(
            gets("sp_opt_dx"),
            gets("sp_opt_xyz"),
            gets("sp_opt_ang"),  # filenames
            dx_dict,  # data about the dx-file (header and how to save it)
            getb("sp_correct"),
            getb("sp_remove"),
            getf("maxval"),  # data concerning postprocessing of energy data
            getb("globalopt"),  # is the global optimum desired?
            obmol,
            np_grid,  # data needed to print out xyz-files at the optimum geometries
            transrot_result,  # see above
        )
