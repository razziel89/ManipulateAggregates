import sys
import os
import copy
import re
from multiprocessing import Pool, Event

import numpy as np

from openbabel import doubleArray, OBAggregate, OBForceField
import pybel as p
from manipulate_molecules import read_from_file
from collection.write import print_dx_file
from collection.read import read_dx

#helper variables for priting of progress output
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

#global data to allow easy sharing of data between processes
global data
data = None #initialized to none
global grid

#this allows for easy data sharing between processes without pickling
def transrot_parallel_init(obmol, transgrid, terminating):
    global data
    data = (obmol, transgrid, terminating)

def _double_array(mylist):
    """Create a C array of doubles from a list."""
    c = doubleArray(len(mylist))
    for i,v in enumerate(mylist):
        c[i] = v
    return c

def _double_dist(iterable):
    """
    Returns a list of (array,-array) where array is a
    C-array made from a single element of iterable.
    
    iterable: numpy array of (float,float,float)
        A numpy array of 3D-vectors that are to be converted
        to plain C-arrays.
    """
    return [(_double_array(i),_double_array(-i)) for i in iterable]

def _gen_trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report,reportstring):
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

def trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report=False,reportstring=""):
    return list(_gen_trans_en(obmol,obff,double_grid,maxval,cutoff,vdw_scale,report,reportstring))

def _print_dx_file(prefix,dictionary,values,comment):
    """
    Helper function to make wring a DX-file easier and ot reduce the
    number of arguments that has to be passed from function to function.
    """
    filename = prefix+dictionary["filename"]
    counts   = dictionary["counts"]
    org      = dictionary["org"]
    delx     = dictionary["delx"]
    dely     = dictionary["dely"]
    delz     = dictionary["delz"]
    gzipped  = dictionary["gzipped"]
    print_dx_file(filename,counts,org,delx,dely,delz,values,comment=comment,gzipped=gzipped)

def _transrot_en_process(args):
    global data

    defaultobmol, transgrid, terminating = data

    try:
        if not terminating.is_set():

            a1, a2, a3, ffname, report, maxval, dx_dict, correct, savetemplate, \
                templateprefix, anglecount, count, save_noopt, save_opt, optsteps, cutoff, vdw_scale, \
                oldfile = args

            angle_string=str(a1)+","+str(a2)+","+str(a3)
            angle_comment="angles=("+angle_string+")"

            compute = True
            if oldfile is not None:
                old                    = read_dx(oldfile,grid=False,data=True,silent=True,comments=True,gzipped=dx_dict["gzipped"])
                old_a1,old_a2,old_a3   = list(map(float,re.split(r',|\(|\)',old["comments"][0])[1:4]))
                if not (a1,a2,a3) == (old_a1,old_a2,old_a3):
                    print >>sys.stderr,"WARNING: old dx-file %s treated %s with index %d. This is also my index but I treat %s. Recomputing."%(oldfile,old["comments"][0],anglecount,angle_comment)
                    compute = True
                else:
                    energies = old['data'].tolist()
                    compute = False
                    del old
                if not len(transgrid) == len(energies):
                    print >>sys.stderr,"WARNING: old dx-file %s contains %d entries but the spatial grid is supposed to have %d entries. Recomputing."%(oldfile,len(energies),len(transgrid))
                    compute = True

            if compute or savetemplate:
                obmol = OBAggregate(defaultobmol)
                obff  = OBForceField.FindForceField(ffname)
                rotfunc=obmol.RotatePart
                rotfunc(0,1,a1)
                rotfunc(0,2,a2)
                rotfunc(0,3,a3)

            if compute:
                energies = trans_en(obmol,obff,transgrid,maxval*1.2,cutoff,vdw_scale,report=report)

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
                    _print_dx_file(str(anglecount)+"_",dx_dict,tempenergies,angle_comment)
            
                if correct or dx_dict["save_dx"]:
                    del tempenergies
            
            if savetemplate:
                minindex = energies.index(min(energies))
                template = grid[minindex][0]
            
                obmol.TranslatePart(0,template)
                if obmol.IsGoodVDW():
                    if save_noopt:
                        filename=templateprefix+str(anglecount)+"_"+angle_string+".xyz"
                        p.Molecule(obmol).write("xyz",filename,overwrite=True)
                    if save_opt:
                        filename=templateprefix+"opt_"+str(anglecount)+"_"+angle_string+".xyz"
                        p_tempmol=p.Molecule(obmol)
                        p_tempmol.localopt(forcefield=ffname,steps=optsteps)
                        p_tempmol.write("xyz",filename,overwrite=True)

            #returning the molecule to its original state is not necessary since every worker process
            #creates its own instance and leaves the original one as is

    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

    return anglecount,(a1,a2,a3),energies,minindex

def transrot_en(obmol,              ffname,
                transgrid,          rotgrid,
                maxval,             dx_dict,        correct,
                cutoff,             vdw_scale,
                report=0,       
                reportcount=1,      reportmax=None,
                savetemplate=True,  templateprefix="template_",
                save_noopt=True,    save_opt=True,     optsteps=500,
                olddxfiles={}
                ):

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
    
    pool = Pool(nr_threads, initializer=transrot_parallel_init, initargs=(obmol, transgrid, terminating))   #NODEBUG
    #global data                            #DEBUG
    #data = (obmol, transgrid, terminating) #DEBUG

    nr_angles = len(rotgrid)
    nr_points = len(transgrid)

    if len(olddxfiles) == 0:
        keyfunc = lambda c,d: d
    else:
        keyfunc = olddxfiles.get

    args=[[a1, a2, a3, 
        ffname, herereport, maxval, dx_dict, correct, 
        savetemplate, templateprefix, anglecount, count, 
        save_noopt, save_opt, optsteps, cutoff, vdw_scale, keyfunc(anglecount,None)] 
        for (a1,a2,a3),anglecount,count in zip(rotgrid,xrange(reportcount,nr_angles+reportcount),xrange(nr_angles))]

    #pre-declare variables
    #the optimum energies
    opt_energies = np.ones((nr_points,),dtype=float)*maxval
    #the angles of the geometries corresponding to the optimum energies
    #360 is the hard-coded default
    opt_angles   = np.ones((nr_points,3),dtype=float)*360.0
    #an element is True if at the corresponding spatial point at least
    #one energy smaller than maxval has been found
    opt_present  = np.zeros((nr_points,),dtype=bool)
    #save the optimum index for each angular arrangement
    opt_angindex = np.zeros((nr_angles,),dtype=int)
    #save the index of the angular arrangement that is the current optimum at the spatial point
    opt_spindex  = np.zeros((nr_points,),dtype=int)
    #an helper array for the comparison
    np_compare   = np.zeros((nr_points,),dtype=bool)
    anglecount=reportcount
    try:
        reportstring = "0/"+str(reportmax)
        print reportstring
        #The structure of temp is: anglecount,(a1,a2,a3),energies,minindex
        #The function _transrot_en_process is guarantueed to return values smaller than maxval only if
        #an actual evaluation using a force field has been performed
        for temp in pool.imap_unordered(_transrot_en_process, args):    #NODEBUG
        #for temp in map(_transrot_en_process, args):                   #DEBUG
            #transform energies to numpy array
            opt_temp = np.array(temp[2])
            #save the optimum index of this angular arrangement for later use
            opt_angindex[temp[0]-reportcount]=temp[3]
            #get all positions where the new energies are smaller than the last optimum
            np.less(opt_temp, opt_energies, out=np_compare)
            #asign those new energies at every point where the comparison was true
            opt_energies[np_compare] = opt_temp[np_compare]
            #asign the angles at every such point
            opt_angles[np_compare] = np.array(temp[1])
            #find out where at least one such assignment has been performed
            opt_present[np_compare] = True
            #which angular arrangement is the current optimum at this spatial point (index)
            opt_spindex[np_compare] = temp[0]-reportcount
            #result.append(temp)
            reportstring = str(anglecount)+"/"+str(reportmax)
            if report!=0:
                print reportstring
            anglecount+=1
        pool.close()
    except KeyboardInterrupt:
        print >>sys.stderr,"Caught keyboard interrupt."
        pool.terminate()    #NODEBUG
        print >>sys.stderr,"Terminating main routine prematurely."
    finally:
        pool.join() #NODEBUG
        #pass       #DEBUG
    return opt_energies,opt_angles,opt_spindex,opt_present,opt_angindex

def _no_none_string(string):
    if string=="None":
        return False
    else:
        return True

def sp_opt(dx, xyz, ang, dx_dict, correct, remove, maxval, globalopt, obmol, grid, transrot_result):
    """
    This function selects optimum geometries and energies for all
    spatial coordinates. It can also sort out such geometries that clashed
    or where the molecules were too far apart from each other.
    """
    dx_bool  = _no_none_string(dx)
    xyz_bool = _no_none_string(xyz)
    ang_bool = _no_none_string(ang)

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
        writefile=p.Outputfile("xyz",xyz,overwrite=True)

        filename=xyz
        tempmol=OBAggregate(obmol)
        pybeltempmol=p.Molecule(tempmol)
        rotfunc=tempmol.RotatePart
        transfunc=tempmol.TranslatePart
        commentfunc=tempmol.SetTitle

        tempgrid = np.copy(grid)
        tempgrid[np.logical_not(opt_present)] = np.array([0.0,0.0,0.0])

        if remove:
            iterable = zip(opt_angles[opt_present],opt_energies[opt_present],tempgrid[opt_present])
        else:
            iterable = zip(opt_angles,opt_energies,tempgrid)

        vec = _double_array([0.0,0.0,0.0])
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
            actualmax = np.amax(opt_energies[opt_present])
            values = np.ones(opt_energies.shape,dtype=float)*actualmax
            values[opt_present] = opt_energies[opt_present]
        else:
            values = opt_energies
        _print_dx_file("",dx_dict,values,"Optimum energies for all spatial points.")

    if globalopt:
        minindex = np.argmin(opt_energies)
        minvalue = opt_energies[minindex]
        mina1,mina2,mina3 = opt_angles[minindex]
        minvec = _double_array(grid[minindex])
        tempmol = OBAggregate(obmol)
        tempmol.RotatePart(0,1,mina1)
        tempmol.RotatePart(0,2,mina2)
        tempmol.RotatePart(0,3,mina3)
        tempmol.TranslatePart(0,minvec)
        tempmol.SetTitle("Energy: "+str(minvalue))
        p.Molecule(tempmol).write("xyz","globalopt.xyz",overwrite=True)
        del tempmol

def general_grid(org,countspos,countsneg,dist,postprocessfunc=None,resetval=False):
    """
    Return a 3D-grid.

    org: the origin of the grid
    countspos: how many points shall be taken in positive x,y and z directions
    countsneg: same as countspos but for negative directions
    dist: the distance for x,y and z directions
    postprocessfunc: if not None, this function is applied to the grid prior
                     to returning it
    resetval: if Ture, return (grid,reset_val) instead of just grid. Here, 
              reset_val is the last element of the original grid that has been
              removed from teh grid an can be used to return the vector to org.
              Meaning adding all elements in grid to org leaves the vector at
              org if resetval==False and at org-reset_val if resetval==True.
    """
    #create helper vectors
    start = org - countsneg*dist
    end   = org + countspos*dist 
    #create grid for rotation of the molecule
    space     = [np.linspace(s,e,num=cp+cn+1,dtype=float) 
                 for s,e,cp,cn 
                 in zip(start,end,countspos,countsneg)
                ]
    a1,a2,a3  = np.array(np.meshgrid(*space,indexing="ij"))
    a1.shape  = (-1,1)
    a2.shape  = (-1,1)
    a3.shape  = (-1,1)
    grid      = np.concatenate((a1,a2,a3),axis=1)
    if postprocessfunc is not None:
        grid = postprocessfunc(grid)
    if resetval:
        reset = -grid[len(grid)-1]
        return grid,reset
    else:
        return grid

def _prepare_molecules(mol1,mol2,aligned_suffix="",save_aligned=False,align=True):
    """
    First, align mol1 and mol2 with their centers to "org"
    and their third and second main axes with [1,0,0] and [0,1,0],
    respectively.
    Then, append mol1 to mol2 and initialize the forcefield
    of the combined OBAggregate-object.
    The OBAggregate object is returned.
    """
    nr_scan_mols=mol2.mol.GetNrMolecules()
    if align:
        #align the molecule's longest axis with the x-axis and the second longest axis with the y-direction
        #and center the molecule to the origin
        mol1.align([0,0,0],[1,0,0],[0,1,0])
        mol2.align([0,0,0],[1,0,0],[0,1,0])
        if save_aligned:
            mol1.write(mol1.fileinfo['name']+aligned_suffix)
            mol2.write(mol2.fileinfo['name']+aligned_suffix)
    #append the molecule
    mol2.append(mol1)
    del(mol1)
    #since python needs quite some time to access an objects member, saving a member saves time
    obmol = mol2.mol
    #configure the aggregate for using tags in order to be able to move multiple
    #molcules in the aggregate with one command
    obmol.EnableTags()
    tag = obmol.CreateTag(nr_scan_mols)
    for i in xrange(nr_scan_mols):
        obmol.AddToTag(i,tag)
    return obmol

def get_old_dxfiles(olddirs,suffix):
    olddxfiles = {}
    dxregex = re.compile("^[1-9][0-9]*_%s$"%(suffix))
    for d in olddirs:
        oldlength = len(olddxfiles)
        if os.path.isdir(d):
            olddxfiles.update({int(f.split("_")[0]):d+os.sep+f for f in os.listdir(d) if re.match(dxregex,f)})
            if len(olddxfiles) == oldlength:
                print >>sys.stderr,"WARNING: directory supposed to contain dx files from previous runs %s does not contain anything matching ^[1-9][0-9]*_%s$ . Skipping."%(d,gets("suffix"))
    return olddxfiles

def scan_main(parser):
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
    for check_option in ["columns","optsteps","progress"]:
        geti(check_option)
    #check whether some options conflict
    if gets("volumetric_data").startswith("from_scan,") and "minimasearch" in gets("jobtype").split(",") and not getb("save_dx"):
        print >>sys.stderr,"WARNING: a subsequent minimasearch tries to get its dx-files from this scan but"
        print >>sys.stderr,"         you requested not to save dx-files. This is probably an error (but not so if"
        print >>sys.stderr,"         you requested those dx-files to be used from a different directory) so please check."

    #value for progress reports
    if geti("progress") not in [0,1,2]:
        raise ValueError('Wrong value for parameter "progress" given. Must be 0,1 or 2.')

    #populate all variables with the given values
    #read in the two molecules/aggregates from the given files
    mol1 = read_from_file(gets("geometry1"),ff=None)
    mol2 = read_from_file(gets("geometry2"),ff=None)
    #spatial grid: check gridtype and set-up grid
    option = "sp_gridtype"
    if gets("sp_gridtype") == "full":
        #these are only the counts in one direction
        np_counts = np.array(map(int,gets("countsxyz").split(",")))
        #example: 0.35,0.5,0.5
        np_del    = np.array(map(float,gets("distxyz").split(",")))
        np_org    = np.array([0,0,0])
        if do_calculate:
            np_grid   = general_grid(np_org,np_counts,np_counts,np_del)
            dx_dict = {"filename": gets("suffix"), "counts": list(2*np_counts+1), "org": list(np_grid[0]),
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
        countsposmain = np.array(map(int,gets("countspos").split(",")))
        countsnegmain = np.array(map(int,gets("countsneg").split(",")))
        distmain      = np.array(map(float,gets("dist").split(",")))
        if do_calculate:
            np_rot        = general_grid(np.array([0.0,0.0,0.0]),countsposmain,countsnegmain,distmain)
    else:
        raise ValueError("Wrong value for config value ang_gridtype.")

    if not do_calculate:
        return

    #align the two molecules and append one to the other
    #after this, mol1 and mol2 can no longer be used
    obmol = _prepare_molecules(mol1,mol2,gets("aligned_suffix"),save_aligned=getb("save_aligned"),align=getb("prealign"))

    #convert the grid to C data types
    grid      = _double_dist(np_grid)

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
    transrot_result = transrot_en(
                obmol,                       gets("forcefield").lower(),
                grid,                        np_rot,
                getf("maxval"),              dx_dict,   getb("correct"),
                getf("cutoff"),              getf("vdw_scale"),
                report=geti("progress"),     reportmax=len(np_rot),
                save_noopt=getb("save_noopt"),
                save_opt=getb("save_opt"),   optsteps=geti("optsteps"),
                olddxfiles=olddxfiles
                )

    del grid #the grid in C data types is no longer needed since the scan has already been performed

    #Evaluate transrot_result to find the angular optimum for every
    #spatial grid point, if so desired
    if getb("sp_opt"):

        dx_dict["filename"]=gets("sp_opt_dx")
        dx_dict["save_dx"]=getb("sp_opt")
        
        sp_opt(
               gets("sp_opt_dx"),   gets("sp_opt_xyz"), gets("sp_opt_ang"), #filenames
               dx_dict, #data about the dx-file (header and how to save it)
               getb("sp_correct"), getb("sp_remove"),  getf("maxval"), #data concerning postprocessing of energy data
               getb("globalopt"), #is the global optimum desired?
               obmol, np_grid, #data needed to print out xyz-files at the optimum geometries
               transrot_result #see above
               )
