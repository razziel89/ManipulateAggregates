import sys
import os
from multiprocessing import Pool, Event

import numpy as np

from openbabel import doubleArray, OBAggregate, OBForceField
import pybel as p
from manipulate_molecules import *
from collection.write import print_dx_file

#example command-line:
#python energyscan.py MKP48.xyz MKP48.xyz mmff94 3 10,10,10 2.5,2.5,2.5 aligned.xyz mmff94.dx 100000 6,6,6 5,5,5 30.0,30.0,30.0 True True
#this will use the geometry in MKP48.xyz, use the force field mmff94
#with a spacing of 2.5 Angstroms in every direction, in total 10 points
#in every direction will be used
#the aligned geometry will be saved to aligned.xyz
#and the scanning result will be savec to mmff94.dx
#for internal computation, if the two molecules clash a value of 100000 will be used
#which will be set to the actual maximum value (the second to last True)
#and progress will be printed every 1000 points (the last True)

#helper variables for priting of progress output
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

#global data to allow easy sharing of data between processes
global data
data = None #initialized to none

#this allows for easy data sharing between processes without pickling
def parallel_init(obmol, transgrid, terminating):
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

def _write_obmol(mol,filename):
     p.Molecule(mol).write("xyz",filename,overwrite=True)

def _gen_trans_en(obmol,obff,double_grid,maxval,report,reportstring):
    if report:
        print ERASE_LINE+"   %s  %.2f%%"%(reportstring,0.0/len(grid))+CURSOR_UP_ONE
    count = 0
    transfunc  = obmol.TranslatePart
    vdwfunc    = obmol.IsGoodVDW
    setupfunc  = obff.Setup
    energyfunc = obff.Energy
    for newvec,retvec in double_grid:
        transfunc(0,newvec)
        if vdwfunc():
            setupfunc(obmol)
            yield energyfunc()
        else:
            yield maxval
        count+=1
        if report and count%1000 == 0:
            print ERASE_LINE+"   %s  %.2f%%"%(reportstring,100.0*count/len(grid))+CURSOR_UP_ONE
        transfunc(0,retvec)
    if report:
        print ERASE_LINE+"   %s  %.2f%%"%(reportstring,100.0)+CURSOR_UP_ONE

def trans_en(obmol,obff,double_grid,maxval,report=False,reportstring=""):
    return list(_gen_trans_en(obmol,obff,double_grid,maxval,report,reportstring))

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

            a1, a2, a3, ffname, report, maxval, dx_dict, correct, savetemplate, templateprefix, anglecount, count = args

            obmol = OBAggregate(defaultobmol)
            obff  = OBForceField.FindForceField(ffname)

            angle_string=str(a1)+","+str(a2)+","+str(a3)
            angle_comment="angles=("+angle_string+")"
            
            rotfunc=obmol.RotatePart
            rotfunc(0,1,a1)
            rotfunc(0,2,a2)
            rotfunc(0,3,a3)

            energies = trans_en(obmol,obff,transgrid,maxval,report)

            if correct:
                try:
                    actualmax = max((e for e in energies if not e==maxval))
                except ValueError:
                    actualmax = maxval
                energies = [actualmax if e==maxval else e for e in energies]
            
            _print_dx_file(str(anglecount)+"_",dx_dict,energies,angle_comment)
            
            if savetemplate:
                minindex = energies.index(min(energies))
                template = grid[minindex][0]
            
                obmol.TranslatePart(0,template)
                if obmol.IsGoodVDW():
                    _write_obmol(obmol,templateprefix+str(anglecount)+"_"+angle_string+".xyz")
            
            #returning the molecule to its original state is not necessary since every worker process
            #creates its own instance and leaves the original one as is

    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

    return

def transrot_en(obmol,              ffname,
                transgrid,          rotgrid,
                maxval,             dx_dict,        correct,
                report=0,       
                reportcount=1,      reportmax=None,
                savetemplate=True,  templateprefix="template_"
                ):
    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
        os.environ["OMP_NUM_THREADS"]=str(1)
    except KeyError:
        os.environ["OMP_NUM_THREADS"]=str(1)
        nr_threads = int(os.environ["OMP_NUM_THREADS"])

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
    
    pool = Pool(nr_threads, initializer=parallel_init, initargs=(obmol, transgrid, terminating))

    args=[[a1, a2, a3, ffname, herereport, maxval, dx_dict, correct, savetemplate, templateprefix, anglecount, count] for (a1,a2,a3),anglecount,count in zip(rotgrid,xrange(reportcount,len(rotgrid)+reportcount),xrange(len(rotgrid)))]

    anglecount=reportcount
    try:
        reportstring = "0/"+str(reportmax)
        print reportstring
        for temp in pool.imap_unordered(_transrot_en_process, args):
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
                 in  zip(start,end,countspos,countsneg)
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
    mol1.write(mol1.fileinfo['name']+"."+aligned_suffix)
    mol2.write(mol2.fileinfo['name']+"."+aligned_suffix)
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

if __name__ == "__main__":

    #treat command-line parameters 
    ffname       = sys.argv[3]
    if ffname.lower() not in ["uff", "mmff94", "gaff", "ghemical"]:
        raise ValueError('Wrong foce field given. Only "uff", "mmff94", "gaff" and "ghemical" will be accepted.')
    #read in the two molecules/aggregates from the given files
    mol1         = read_from_file(sys.argv[1],ff=None)
    mol2         = read_from_file(sys.argv[2],ff=None)
    #this is the number of coloumns to be printed to the dx-file
    coloumns  = int(sys.argv[4])
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
    if len(sys.argv)>=14:
        if sys.argv[13]=="1":
            correct=True
        else:
            correct=False
    else:
        correct=False

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
    
    obmol                 = _prepare_molecules(mol1,mol2,aligned_suffix)

    np_grid,np_reset_vec  = general_grid(np_org,np_counts,np_counts,np_del,resetval=True)

    grid                  = _double_dist(np_grid)

    np_rot                = general_grid(np.array([0.0,0.0,0.0]),countsposmain,countsposmain,distmain)

    dx_dict = {"filename": dx_file, "counts": list(2*np_counts+1), "org": list(np_grid[0]),
               "delx": [np_del[0],0.0,0.0], "dely": [0.0,np_del[1],0.0], "delz": [0.0,0.0,np_del[2]]}

    transrot_en(obmol,  ffname,
                grid,   np_rot,
                maxval, dx_dict,        correct,
                report=report,          reportmax=len(np_rot),
                )
