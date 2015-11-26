import sys
import os
import copy
from multiprocessing import Pool, Event

import numpy as np

from openbabel import doubleArray, OBAggregate, OBForceField
import pybel as p
from manipulate_molecules import *
from collection.write import print_dx_file

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
    Returns a list of double arrays each of which contains
    the list-difference to the preceeding list as well as
    a double array to reset back to [0,0,0].
    E.g. the input
        [ [ 1, 2, 3],
          [ 4, 5, 6], 
          [ 7, 8, 9] ]
        would return a list of double arrays containing
        [ [ 1, 2, 3],
          [ 3, 3, 3],
          [ 3, 3, 3],
          [-7,-8,-9] ]
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
    print_dx_file(filename,counts,org,delx,dely,delz,values,comment=comment)

def _transrot_en_process(args):
    global data

    defaultobmol, transgrid, terminating = data

    try:
        if not terminating.is_set():

            a1, a2, a3, ffname, report, maxval, dx_dict, correct, savetemplate, \
                templateprefix, anglecount, count, save_noopt, save_opt, optsteps, cutoff, vdw_scale = args

            obmol = OBAggregate(defaultobmol)
            obff  = OBForceField.FindForceField(ffname)

            angle_string=str(a1)+","+str(a2)+","+str(a3)
            angle_comment="angles=("+angle_string+")"
            
            rotfunc=obmol.RotatePart
            rotfunc(0,1,a1)
            rotfunc(0,2,a2)
            rotfunc(0,3,a3)

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
                ):
    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
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
    
    pool = Pool(nr_threads, initializer=transrot_parallel_init, initargs=(obmol, transgrid, terminating))

    nr_angles = len(rotgrid)
    nr_points = len(transgrid)

    args=[[a1, a2, a3, 
        ffname, herereport, maxval, dx_dict, correct, 
        savetemplate, templateprefix, anglecount, count, 
        save_noopt, save_opt, optsteps, cutoff, vdw_scale] 
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
        for temp in pool.imap_unordered(_transrot_en_process, args):
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
        pool.terminate()
        print >>sys.stderr,"Terminating main routine prematurely."
    finally:
        pool.join()
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
        tempmol=op.OBAggregate(obmol)
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
        tempmol = op.OBAggregate(obmol)
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

def _prepare_molecules(mol1,mol2,aligned_suffix):
    """
    First, align mol1 and mol2 with their centers to "org"
    and their third and second main axes with [1,0,0] and [0,1,0],
    respectively.
    Then, append mol1 to mol2 and initialize the forcefield
    of the combined OBAggregate-object.
    The OBAggregate and OBForcefield objects are returned.
    """
    nr_scan_mols=mol2.mol.GetNrMolecules()
    #align the molecule's longest axis with the x-axis and the second longest axis with the y-direction
    #and center the molecule to the origin
    mol1.align([0,0,0],[1,0,0],[0,1,0])
    mol2.align([0,0,0],[1,0,0],[0,1,0])
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

def _bool_parameter(index, default):
    if len(sys.argv)>=index+1:
        if sys.argv[index]=="1":
            result=True
        else:
            result=False
    else:
        result=default
    return result

def print_example():
    """
    Print an example config file to stdout.
    """
    s="""#all lines starting with # are comments and can be removed
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

def newmain():
    global grid
    #default configuration values
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

    from ConfigParser import NoOptionError
    from collection.read import read_config_file as rf
    parser = rf(sys.argv[1],defaults=config)
    gets   = parser.get_str
    geti   = parser.get_int
    getf   = parser.get_float
    getb   = parser.get_boolean
    #do some error checking
    #forcefield
    if gets("forcefield").lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    #value for progress reports
    if geti("progress") not in [0,1,2]:
        raise ValueError('Wrong value for parameter "progress" given. Must be 0,1 or 2.')
    #boolean values
    for check_option in ["save_dx","save_aligned","save_noopt","save_opt","correct","sp_opt"]:
        getb(check_option) #this will throw errors if the value cannot be converted to a boolean
    #remaining float values
    try:
        getf("cutoff")
    except ValueError:
        raise TypeError("Option cutoff must be of type float.")
    try:
        getf("vdw_scale")
    except ValueError:
        raise TypeError("Option vdw_scale must be of type float.")
    #remaining integer values
    try:
        geti("columns")
    except ValueError:
        raise TypeError("Option cutoff must be of type int.")
    #check whether some options conflict
    #NO CONFLICTS KNOWN YET
    #populate all variables with the given values
    try:
        #read in the two molecules/aggregates from the given files
        option="geometry1"
        mol1 = read_from_file(gets(option),ff=None)
        option="geometry2"
        mol2 = read_from_file(gets(option),ff=None)
        #spatial grid: check gridtype and set-up grid
        if gets("sp_gridtype") == "full":
            #these are only the counts in one direction
            option="countsxyz"
            np_counts = np.array(map(int,gets(option).split(",")))
            #example: 0.35,0.5,0.5
            option="distxyz"
            np_del    = np.array(map(float,gets(option).split(",")))
            np_org    = np.array([0,0,0])
            np_grid   = general_grid(np_org,np_counts,np_counts,np_del)
            option="suffix"
            dx_dict = {"filename": gets(option), "counts": list(2*np_counts+1), "org": list(np_grid[0]),
                       "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}
            option="save_dx"
            dx_dict[option]=getb(option)
        else:
            raise ValueError("Wrong value for config value sp_gridtype.")
        #angular grid: check gridtype and set-up grid
        if gets("ang_gridtype") == "full":
            #these are the counts and distances for rotation
            option="countspos"
            countsposmain = np.array(map(int,gets(option).split(",")))
            option="countsneg"
            countsnegmain = np.array(map(int,gets(option).split(",")))
            option="dist"
            distmain      = np.array(map(float,gets(option).split(",")))
            np_rot        = general_grid(np.array([0.0,0.0,0.0]),countsposmain,countsposmain,distmain)
        else:
            raise ValueError("Wrong value for config value ang_gridtype.")
    except NoOptionError:
        raise KeyError("Necessary option missing from config file: "+option)

    #align the two molecules and append one to the other
    #after this, mol1 and mol2 can no longer be used
    obmol = _prepare_molecules(mol1,mol2,gets("aligned_suffix"))

    #convert the grid to C data types
    grid      = _double_dist(np_grid)

    #For every angle, scan the entire spatial grid and save
    #each optimum geometry if desired
    #Will also return a structure making it easy to find the optimum
    #for every spatial point
    transrot_result = transrot_en(
                obmol,                       gets("forcefield"),
                grid,                        np_rot,
                getf("maxval"),              dx_dict,   getb("correct"),
                getf("cutoff"),              getf("vdw_scale"),
                report=geti("progress"),     reportmax=len(np_rot),
                save_noopt=getb("save_noopt"),
                save_opt=getb("save_opt"),   optsteps=geti("optsteps")
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

def oldmain():
    global grid
    print >>sys.stderr, "WARNING: providing options on the command line is deprecated, use config file instead."

    cutoff=25.0 #hard-coded so far
    #treat command-line parameters 
    ffname       = sys.argv[3]
    if ffname.lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    #read in the two molecules/aggregates from the given files
    mol1         = read_from_file(sys.argv[1],ff=None)
    mol2         = read_from_file(sys.argv[2],ff=None)
    #this is the number of columns to be printed to the dx-file
    columns  = int(sys.argv[4])
    #example: 40,20,20
    #these are only the counts in one direction
    np_counts = np.array(map(int,sys.argv[5].split(",")))
    #example: 0.35,0.5,0.5
    np_del    = np.array(map(float,sys.argv[6].split(",")))
    np_org    = np.array([0,0,0])
    
    #some additional parameters
    aligned_suffix = sys.argv[7]
    dx_file        = sys.argv[8]
    maxval         = float(sys.argv[9])
    
    #these are the counts and distances for rotation
    countsposmain = np.array(map(int,sys.argv[10].split(",")))
    countsnegmain = np.array(map(int,sys.argv[11].split(",")))
    distmain      = np.array(map(float,sys.argv[12].split(",")))
    
    #whether or not the output shall be corrected so that
    #geometries with VDW-clashes are set to the maximum energy
    #of all energies for which no clashes occurred
    correct = _bool_parameter(13,False)

    #whether or not progress reports shall be performed
    if len(sys.argv)>=15:
        if sys.argv[14]=="1":
            report=1
        elif sys.argv[14]=="2":
            report=2
        else:
            report=0
    else:
        report=0

    save_noopt = _bool_parameter(15,True)
    save_opt   = _bool_parameter(16,False)

    if len(sys.argv)>=18:
        optsteps=int(sys.argv[17])
    else:
        optsteps = 500

    obmol                 = _prepare_molecules(mol1,mol2,aligned_suffix)

    np_grid,np_reset_vec  = general_grid(np_org,np_counts,np_counts,np_del,resetval=True)
    grid                  = _double_dist(np_grid)

    np_rot                = general_grid(np.array([0.0,0.0,0.0]),countsposmain,countsposmain,distmain)

    dx_dict = {"filename": dx_file, "counts": list(2*np_counts+1), "org": list(np_grid[0]),
               "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}

    transrot_en(obmol,                  ffname,
                grid,                   np_rot,
                maxval,                 dx_dict,                correct,
                cutoff,                 -1.0, #this is the vdw scale and hardcoded
                report=report,          reportmax=len(np_rot),
                save_noopt=save_noopt,  save_opt=save_opt,      optsteps=optsteps
                )

if __name__ == "__main__":
    if len(sys.argv)==1:
        print_example()
    elif len(sys.argv)==2:
        newmain()
    elif len(sys.argv)>13:
        oldmain()
    else:
        raise ValueError("Wrong number of parameters given. Launch without parameters to get an example config file.")
