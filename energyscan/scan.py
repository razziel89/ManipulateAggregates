"""Systematically determine aggregate energies.

This subsubmodule is part of ManipulateAggregates.energyscan. It implements the
first step of the 3-step procedure that creates low energy aggregate geometries.

Parallelization is supported for this subsubmodule.

@package ManipulateAggregates.energyscan.scan
"""

#This file is part of ManipulateAggregates.
#
#Copyright (C) 2016 by Torsten Sachse
#
#ManipulateAggregates is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ManipulateAggregates is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import copy
import re
import errno
from multiprocessing import Pool, Event

import logging
logger = logging.getLogger(__name__)
try:
    import numpy as numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from openbabel import doubleArray, OBAggregate, OBForceField
except ImportError:
    logger.warning("Could not import doubleArray, OBAggregate or OBForceField from openbabel")
try:
    import pybel as pybel
except ImportError:
    logger.warning("Could not import pybel")
try:
    from ..aggregate import read_from_file
except ImportError:
    logger.warning("Could not import ..aggregate.read_from_file")
try:
    from ..energyscan.ansilliary import CURSOR_UP_ONE, ERASE_LINE, init_hashing, double_dist, double_array
except ImportError:
    logger.warning("Could not import CURSOR_UP_ONE, ERASE_LINE, init_hashing, double_dist or double_array from ..energyscan.ansilliary")
try:
    from ..energyscan.ansilliary import print_dx_file, general_grid, prepare_molecules, get_old_dxfiles, no_none_string
except ImportError:
    logger.warning("Could not import print_dx_file, general_grid, prepare_molecules, get_old_dxfiles or no_none_string from ..energyscan.ansilliary")
try:
    from ..collection.read import read_dx
except ImportError:
    logger.warning("Could not import ..collection.read.read_dx")
try:
    from ..collection import hashIO
except ImportError:
    logger.warning("Could not import ..collection.read.hashIO")

## \cond
#global data to allow easy sharing of data between processes
global data_s
data_s = None #initialized to none
global grid
## \endcond

#this allows for easy data sharing between processes without pickling
def _transrot_parallel_init(obmol, transgrid, terminating):
    """Allow easy data sharing between processes without pickling.

    Args:
        obmol: (OpenBabel OBMol object): molecule to be treated
        transgrid: (list of C-arrays of type double) spatial grid for the scan
        terminating: (return value of multiprocessing.Event()) specifies 
            whether or not a parallel computation has been terminated
            prematurely or not
    """
    global data_s
    data_s = (obmol, transgrid, terminating)

def _gen_trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report,reportstring):
    """Internal function, undocumented."""
    if report:
        print ERASE_LINE+"   %s  %.2f%%"%(reportstring,0.0/len(grid))+CURSOR_UP_ONE
    count = 0
    transfunc  = obmol.TranslatePart
    #vdwfunc    = obmol.IsGoodVDW
    vdwfunc    = obmol.MinVDWDist
    setupfunc  = obff.Setup
    energyfunc = obff.Energy
    for newvec,retvec in double_grid:
        transfunc(0,newvec)
        dist = vdwfunc(True,vdw_scale) #True will cause openbabel to interrupt whenever a vdW clash was found
        if dist>0.0 and dist<cutoff:
        #if vdwfunc():
            setupfunc(obmol)
            yield energyfunc(False) #this will cause openbabel to not evaluate gradients
        else:
            yield maxval
        count+=1
        if report and count%1000 == 0:
            print ERASE_LINE+"   %s  %.2f%%"%(reportstring,100.0*count/len(grid))+CURSOR_UP_ONE
        transfunc(0,retvec)
    if report:
        print ERASE_LINE+"   %s  %.2f%%"%(reportstring,100.0)+CURSOR_UP_ONE

def _trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report=False,reportstring=""):
    """Internal function, undocumented."""
    return list(_gen_trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report,reportstring))

def _transrot_en_process(args):
    """Each worker process executes this function.

    Args:
        args: (list) arguments to be passed to the worker processes (via
            pickling)
    """
    global data_s

    defaultobmol, transgrid, terminating = data_s

    try:
        if not terminating.is_set():

            a1, a2, a3, ffname, report, maxval, dx_dict, correct, savetemplate, \
                templateprefix, anglecount, count, save_noopt, save_opt, optsteps, cutoff, vdw_scale, \
                oldfile = args

            angle_string=str(a1)+","+str(a2)+","+str(a3)
            angle_comment="angles=("+angle_string+")"

            if oldfile is not None:
                compute = False
                try:
                    old = read_dx(oldfile,grid=False,data=True,silent=True,comments=True,gzipped=dx_dict["gzipped"])
                except ValueError as e:
                    print >>sys.stderr,"Error when reading in old dx-file %s, recomputing. Error was:"%(oldfile),e
                    compute = True
                if not compute:
                    old_a1,old_a2,old_a3   = list(map(float,re.split(r',|\(|\)',old["comments"][0])[1:4]))
                    if not (a1,a2,a3) == (old_a1,old_a2,old_a3):
                        print >>sys.stderr,"WARNING: old dx-file %s treated %s with index %d. This is also my index but I treat %s. Recomputing."%(oldfile,old["comments"][0],anglecount,angle_comment)
                        compute = True
                    else:
                        energies = old['data'].tolist()
                        del old
                if not compute:
                    if not len(transgrid) == len(energies):
                        print >>sys.stderr,"WARNING: old dx-file %s contains %d entries but the spatial grid is supposed to have %d entries. Recomputing."%(oldfile,len(energies),len(transgrid))
                        compute = True
            else:
                compute = True

            if compute or savetemplate:
                obmol = OBAggregate(defaultobmol)
                obff  = OBForceField.FindForceField(ffname)
                rotfunc=obmol.RotatePart
                rotfunc(0,1,a1)
                rotfunc(0,2,a2)
                rotfunc(0,3,a3)

            if compute:
                energies = _trans_en(obmol,obff,transgrid,maxval*1.2,cutoff,vdw_scale,report=report)

                if correct or dx_dict["save_dx"]:
                    #create a copy which can then be changed and possibly saved
                    tempenergies = copy.copy(energies)

                if correct:
                    try:
                        actualmax = max((e for e in tempenergies if not e>=maxval))
                    except ValueError:
                        actualmax = maxval
                    tempenergies = [actualmax if e>=maxval else e for e in tempenergies]
            
                if dx_dict["save_dx"]:
                    print_dx_file(str(anglecount)+"_",True,dx_dict,tempenergies,angle_comment)
            
                if correct or dx_dict["save_dx"]:
                    del tempenergies
            
            if savetemplate:
                minindex = energies.index(min(energies))
                template = grid[minindex][0]
            
                obmol.TranslatePart(0,template)
                if obmol.IsGoodVDW(vdw_scale):
                    if save_noopt:
                        filename=templateprefix+str(anglecount)+"_"+angle_string+".xyz"
                        pybel.Molecule(obmol).write("xyz",filename,overwrite=True)
                    if save_opt:
                        filename=templateprefix+"opt_"+str(anglecount)+"_"+angle_string+".xyz"
                        p_tempmol=pybel.Molecule(obmol)
                        p_tempmol.localopt(forcefield=ffname,steps=optsteps)
                        p_tempmol.write("xyz",filename,overwrite=True)

            #returning the molecule to its original state is not necessary since every worker process
            #creates its own instance and leaves the original one as is

    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

    return anglecount,(a1,a2,a3),energies,minindex

def _transrot_en(obmol,             ffname,
                transgrid,          rotgrid,
                maxval,             dx_dict,        correct,
                cutoff,             vdw_scale,
                report=0,       
                reportcount=1,      reportmax=None,
                savetemplate=True,  templateprefix="template_",
                save_noopt=True,    save_opt=True,     optsteps=500,
                olddxfiles={},      partition=(1,1)
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
            print >>sys.stderr, "You requested detailed progress reports but that is only supported for single threaded"
            print >>sys.stderr, "calculations (you chose "+str(nr_threads)+" threads). Will switch to semi-detailed reports."

    #how to properly handle keyboard interrupts when multi processing has been taken from:
    #http://stackoverflow.com/questions/14579474/multiprocessing-pool-spawning-new-childern-after-terminate-on-linux-python2-7
    terminating = Event()
    
    pool = Pool(nr_threads, initializer=_transrot_parallel_init, initargs=(obmol, transgrid, terminating))   #NODEBUG
    #global data_s                            #DEBUG
    #data_s = (obmol, transgrid, terminating) #DEBUG

    nr_angles = len(rotgrid)
    nr_points = len(transgrid)

    if len(olddxfiles) == 0:
        keyfunc = lambda c,d: d
    else:
        keyfunc = olddxfiles.get

    if partition[1] > nr_angles:
        raise ValueError("Number of partitions %d cannot be greater than the number of angles %d"%(partition[1],nr_angles))

    if not partition == (1,1):
        print "...this is a partitioned calculation (partition %d out of %d)..."%partition

        mypartition            = partition[0]
        angles_per_partition   = nr_angles/(partition[1])
        additional_angles      = nr_angles%(partition[1])
        partition_start_end    = [None]*(partition[1])
        partition_start_end[0] = (0, angles_per_partition + (1 if additional_angles>0 else 0))
        additional_angles      -= 1
        for p in xrange(1,partition[1]):
            partition_start_end[p] = (
                    partition_start_end[p-1][1] ,
                    partition_start_end[p-1][1] + angles_per_partition + (1 if additional_angles>0 else 0)
                    )
            additional_angles      -= 1

        print "...partitions are:"
        for (ps,pe),c in zip(partition_start_end,range(1,partition[1]+1)):
            print "          %2d: START: %6d - END: %6d"%(c,ps+1,pe)+("   <-- that's me" if c==mypartition else "")
        mypartition = (partition_start_end[mypartition-1][0],partition_start_end[mypartition-1][1])
        reportmax = mypartition[1] - mypartition[0]
    else:
        print "...this is no partitioned calculation..."
        partition_start_end = [(0,nr_angles)]
        mypartition         = (0,nr_angles)

    args=[[a1, a2, a3, 
        ffname, herereport, maxval, dx_dict, correct, 
        savetemplate, templateprefix, anglecount, count, 
        save_noopt, save_opt, optsteps, cutoff, vdw_scale, keyfunc(anglecount,None)] 
        for (a1,a2,a3),anglecount,count in zip(rotgrid,xrange(reportcount,nr_angles+reportcount),xrange(nr_angles))
        ][mypartition[0]:mypartition[1]]

    #pre-declare variables
    #the optimum energies
    opt_energies = numpy.ones((nr_points,),dtype=float)*maxval
    #the angles of the geometries corresponding to the optimum energies
    #360 is the hard-coded default
    opt_angles   = numpy.ones((nr_points,3),dtype=float)*360.0
    #an element is True if at the corresponding spatial point at least
    #one energy smaller than maxval has been found
    opt_present  = numpy.zeros((nr_points,),dtype=bool)
    #save the optimum index for each angular arrangement
    opt_angindex = numpy.zeros((nr_angles,),dtype=int)
    #save the index of the angular arrangement that is the current optimum at the spatial point
    opt_spindex  = numpy.zeros((nr_points,),dtype=int)
    #an helper array for the comparison
    np_compare   = numpy.zeros((nr_points,),dtype=bool)
    anglecount   = reportcount
    try:
        if mypartition[0] > 0:
            reportstring = "START --> skipping %d"%(mypartition[0])
        else:
            reportstring = "START"
        print reportstring
        #The structure of temp is: anglecount,(a1,a2,a3),energies,minindex
        #The function _transrot_en_process is guarantueed to return values smaller than maxval only if
        #an actual evaluation using a force field has been performed
        for temp in pool.imap_unordered(_transrot_en_process, args):    #NODEBUG
        #for arg in args:                                               #DEBUG
            #temp = _transrot_en_process(arg)                           #DEBUG
            #transform energies to numpy array
            opt_temp = numpy.array(temp[2])
            #save the optimum index of this angular arrangement for later use
            opt_angindex[temp[0]-reportcount]=temp[3]
            #get all positions where the new energies are smaller than the last optimum
            numpy.less(opt_temp, opt_energies, out=np_compare)
            #asign those new energies at every point where the comparison was true
            opt_energies[np_compare] = opt_temp[np_compare]
            #asign the angles at every such point
            opt_angles[np_compare] = numpy.array(temp[1])
            #find out where at least one such assignment has been performed
            opt_present[np_compare] = True
            #which angular arrangement is the current optimum at this spatial point (index)
            opt_spindex[np_compare] = temp[0]-reportcount
            #result.append(temp)
            reportstring = "%d/%d == #%d"%(anglecount,reportmax,mypartition[0]+anglecount)
            if report!=0:
                print reportstring
            anglecount+=1
        if nr_angles-mypartition[0]-reportmax > 0:
            reportstring = "END --> skipping %d"%(nr_angles-mypartition[0]-reportmax)
        else:
            reportstring = "END"
        print reportstring
        pool.close()    #NODEBUG
        pool.join()     #NODEBUG
    except KeyboardInterrupt as e:
        print >>sys.stderr,"Caught keyboard interrupt."
        pool.terminate()    #NODEBUG
        pool.join()         #NODEBUG
        print >>sys.stderr,"Terminating main routine prematurely."
        raise e

    return opt_energies,opt_angles,opt_spindex,opt_present,opt_angindex


def _sp_opt(dx, xyz, ang, dx_dict, correct, remove, maxval, globalopt, obmol, grid, transrot_result):
    """
    This function selects optimum geometries and energies for all
    spatial coordinates. It can also sort out such geometries that clashed
    or where the molecules were too far apart from each other.
    """
    dx_bool  = no_none_string(dx)
    xyz_bool = no_none_string(xyz)
    ang_bool = no_none_string(ang)

    opt_energies,opt_angles,opt_spindex,opt_present,opt_angindex = transrot_result

    if ang_bool:
        import csv
        #how to write csv files taken from https://docs.python.org/2/library/csv.html
        with open(ang, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',',
                                    quotechar  ='|', quoting=csv.QUOTE_MINIMAL)
            if remove:
                iterable = zip(opt_angles[opt_present],grid[opt_present])
            else:
                iterable = zip(opt_angles,grid)

            for (a1,a2,a3),(x,y,z) in iterable:
                spamwriter.writerow([a1,a2,a3,x,y,z])

    if xyz_bool:
        writefile=pybel.Outputfile("xyz",xyz,overwrite=True)

        filename=xyz
        tempmol=OBAggregate(obmol)
        pybeltempmol=pybel.Molecule(tempmol)
        rotfunc=tempmol.RotatePart
        transfunc=tempmol.TranslatePart
        commentfunc=tempmol.SetTitle

        tempgrid = numpy.copy(grid)
        tempgrid[numpy.logical_not(opt_present)] = numpy.array([0.0,0.0,0.0])

        if remove:
            iterable = zip(opt_angles[opt_present],opt_energies[opt_present],tempgrid[opt_present])
        else:
            iterable = zip(opt_angles,opt_energies,tempgrid)

        vec = double_array([0.0,0.0,0.0])
        for (a1,a2,a3),e,(x,y,z) in iterable:
                commentfunc("Energy: "+str(e))
                rotfunc(0,1,a1)
                rotfunc(0,2,a2)
                rotfunc(0,3,a3)
                vec[0]=x
                vec[1]=y
                vec[2]=z
                tempmol.TranslatePart(0,vec)
                writefile.write(pybeltempmol)
                vec[0]=-x
                vec[1]=-y
                vec[2]=-z
                tempmol.TranslatePart(0,vec)
                rotfunc(0,3,-a3)
                rotfunc(0,2,-a2)
                rotfunc(0,1,-a1)

        del tempmol
        writefile.close()

    if dx_bool:
        if correct:
            actualmax = numpy.amax(opt_energies[opt_present])
            values = numpy.ones(opt_energies.shape,dtype=float)*actualmax
            values[opt_present] = opt_energies[opt_present]
        else:
            values = opt_energies
        print_dx_file("",False,dx_dict,values,"Optimum energies for all spatial points.")

    if globalopt:
        minindex = numpy.argmin(opt_energies)
        minvalue = opt_energies[minindex]
        mina1,mina2,mina3 = opt_angles[minindex]
        minvec = double_array(grid[minindex])
        tempmol = OBAggregate(obmol)
        tempmol.RotatePart(0,1,mina1)
        tempmol.RotatePart(0,2,mina2)
        tempmol.RotatePart(0,3,mina3)
        tempmol.TranslatePart(0,minvec)
        tempmol.SetTitle("Energy: "+str(minvalue))
        pybel.Molecule(tempmol).write("xyz","globalopt.xyz",overwrite=True)
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
    gets   = parser.get_str
    geti   = parser.get_int
    getf   = parser.get_float
    getb   = parser.get_boolean
    do_calculate = not(getb("config_check"))
    #do some error checking
    #forcefield
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    temp_ff = OBForceField.FindType(gets("forcefield").lower())
    if temp_ff is None:
        raise ValueError("Somehow there was an error loading the forcefield %s (although it should be known to OpenBabel)."%(gets("forcefield").lower()))
    del temp_ff
    #boolean values
    for check_option in ["save_dx","save_aligned","save_noopt","save_opt","correct","sp_opt","sp_correct",
            "sp_remove","globalopt","prealign","gzipped"]:
        getb(check_option)
    #remaining float values
    for check_option in ["cutoff","vdw_scale","maxval","cutoff","vdw_scale"]:
        getf(check_option)
    #remaining integer values
    for check_option in ["columns","optsteps","progress","hashwidth","hashdepth"]:
        geti(check_option)
    #check whether some options conflict
    if gets("volumetric_data").startswith("from_scan,") and "minimasearch" in gets("jobtype").split(",") and not getb("save_dx"):
        print >>sys.stderr,"WARNING: a subsequent minimasearch tries to get its dx-files from this scan but"
        print >>sys.stderr,"         you requested not to save dx-files. This is probably an error (but not so if"
        print >>sys.stderr,"         you requested those dx-files to be used from a different directory) so please check."

    #initialize directory name hashing
    init_hashing(geti("hashdepth"),geti("hashwidth"),gets("hashalg"))

    #value for progress reports
    if geti("progress") not in [0,1,2]:
        raise ValueError('Wrong value for parameter "progress" given. Must be 0,1 or 2.')

    #populate all variables with the given values
    #read in the two molecules/aggregates from the given files
    mol1 = read_from_file(gets("geometry1"),ff=None)
    mol2 = read_from_file(gets("geometry2"),ff=None)
    #spatial grid: check gridtype and set-up grid
    if gets("sp_gridtype") == "full":
        #these are only the counts in one direction
        np_counts = numpy.array(map(int,gets("countsxyz").split(",")))
        #example: 0.35,0.5,0.5
        np_del    = numpy.array(map(float,gets("distxyz").split(",")))
        np_org    = numpy.array([0,0,0])
        if do_calculate:
            np_grid   = general_grid(np_org,np_counts,np_counts,np_del)
            dx_dict = {"filename": gets("suffix"), "counts": list(2*np_counts+1), "org": list(np_grid[0]),
                       "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}
            dx_dict["save_dx"]=getb("save_dx")
            dx_dict["gzipped"]=getb("gzipped")
        else:
            gets("suffix")
            getb("save_dx")
    elif gets("sp_gridtype") == "half":
        #these are only the counts in one direction
        np_counts_pos = numpy.array(map(int,gets("countsxyz").split(",")))
        np_counts_neg = numpy.array(map(int,gets("countsxyz").split(",")))
        halfspace_vec = list(map(int,gets("halfspace").split(",")))
        for i in xrange(3):
            if halfspace_vec[i]<0:
                np_counts_pos[i] = abs(halfspace_vec[i])
            if halfspace_vec[i]>0:
                np_counts_neg[i] = abs(halfspace_vec[i])
        np_del    = numpy.array(map(float,gets("distxyz").split(",")))
        np_org    = numpy.array([0,0,0])
        if do_calculate:
            np_grid   = general_grid(np_org,np_counts_pos,np_counts_neg,np_del)
            dx_dict = {"filename": gets("suffix"), "counts": list(np_counts_pos+np_counts_neg+1), "org": list(np_grid[0]),
                       "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}
            dx_dict["save_dx"]=getb("save_dx")
            dx_dict["gzipped"]=getb("gzipped")
        else:
            gets("suffix")
            getb("save_dx")
    else:
        raise ValueError("Wrong value for config value sp_gridtype.")
    #check whether this gives an error
    restarted = (len(gets("scan_restartdirs")) > 0)
    if restarted:
        olddirs = gets("scan_restartdirs").split(",")
        for d in olddirs:
            if not os.path.isdir(d):
                if do_calculate:
                    print >>sys.stderr,"WARNING: directory supposed to contain dx files from previous runs %s does not exist. Skipping."%(d)
                else:
                    raise ValueError("Directory supposed to contain dx files from previous runs %s does not exist."%(d))
    #angular grid: check gridtype and set-up grid
    if gets("ang_gridtype") == "full":
        #these are the counts and distances for rotation
        countsposmain = numpy.array(map(int,gets("countspos").split(",")))
        countsnegmain = numpy.array(map(int,gets("countsneg").split(",")))
        distmain      = numpy.array(map(float,gets("dist").split(",")))
        if do_calculate:
            np_rot        = general_grid(numpy.array([0.0,0.0,0.0]),countsposmain,countsnegmain,distmain)
    else:
        raise ValueError("Wrong value for config value ang_gridtype.")

    partition = tuple(map(int,gets("partition").split("/")))
    if len(partition) != 2:
        raise ValueError("Format for 'partition' must be I1/I2 with I1 and I2 positive integers and I1<=I2")
    if partition[0]>partition[1] or partition[0]<1 or partition[1]<1:
        raise ValueError("Format for 'partition' must be I1/I2 with I1 and I2 positive integers and I1<=I2")

    if not do_calculate:
        return

    #align the two molecules and append one to the other
    #after this, mol1 and mol2 can no longer be used
    obmol = prepare_molecules(mol1,mol2,gets("aligned_suffix"),save_aligned=getb("save_aligned"),align=getb("prealign"))

    #convert the grid to C data types
    grid      = double_dist(np_grid)

    if restarted:
        print "This is a restarted run (old files are in: %s)"%(gets("scan_restartdirs"))
        olddxfiles = get_old_dxfiles(gets("scan_restartdirs").split(","),gets("suffix"))
        print "Number of already existing dx files: %d"%(len(olddxfiles))
    else:
        olddxfiles = {}

    #For every angle, scan the entire spatial grid and save
    #each optimum geometry if desired
    #Will also return a structure making it easy to find the optimum
    #for every spatial point
    transrot_result = _transrot_en(
                obmol,                       gets("forcefield").lower(),
                grid,                        np_rot,
                getf("maxval"),              dx_dict,   getb("correct"),
                getf("cutoff"),              getf("vdw_scale"),
                report=geti("progress"),     reportmax=len(np_rot),
                save_noopt=getb("save_noopt"),
                save_opt=getb("save_opt"),   optsteps=geti("optsteps"),
                olddxfiles=olddxfiles,       partition=partition
                )

    del grid #the grid in C data types is no longer needed since the scan has already been performed

    #Evaluate transrot_result to find the angular optimum for every
    #spatial grid point, if so desired
    if getb("sp_opt"):

        dx_dict["filename"]=gets("sp_opt_dx")
        dx_dict["save_dx"]=getb("sp_opt")
        
        _sp_opt(
               gets("sp_opt_dx"),   gets("sp_opt_xyz"), gets("sp_opt_ang"), #filenames
               dx_dict, #data about the dx-file (header and how to save it)
               getb("sp_correct"), getb("sp_remove"),  getf("maxval"), #data concerning postprocessing of energy data
               getb("globalopt"), #is the global optimum desired?
               obmol, np_grid, #data needed to print out xyz-files at the optimum geometries
               transrot_result #see above
               )
