from manipulate_molecules import *
import numpy as np
from collection.write import print_dx_file
import sys

#example command-line:
#python energyscan.py MKP48.xyz mmff94 3 10,10,10 2.5,2.5,2.5 aligned.xyz mmff94.dx 100000 True True
#this will use the geometry in MKP48.xyz, use the force field mmff94
#with a spacing of 2.5 Angstroms in every direction, in total 10 points
#in every direction will be used
#the aligned geometry will be saved to aligned.xyz
#and the scanning result will be savec to mmff94.dx
#for internal computation, if the two molecules clash a value of 100000 will be used
#which will be set to the actual maximum value (the second to last True)
#and progress will be printed every 1000 points (the last True)

def gen_double_dist(np_grid):
    oldpos    = np.zeros((3),dtype=float)
    for curpos in np_grid:
        yield double_array(curpos-oldpos)
        oldpos = curpos
    yield double_array(-curpos)

def double_dist(np_grid):
    return list(gen_double_dist(np_grid))

def gen_trans_en(obmol,obff,double_grid,reset_vec,maxval,report=False):
    if report:
        print "%d/%d"%(0,len(grid))
    count = 0
    transfunc = obmol.TranslatePart
    for newvec in double_grid:
        transfunc(0,newvec)
        if obmol.IsGoodVDW():
            obff.Setup(obmol)
            yield obff.Energy()
        else:
            yield maxval
        count+=1
        if report and count%1000 == 0:
            print "%d/%d"%(count,len(grid))
    transfunc(0,reset_vec)
    if report:
        print "%d/%d"%(len(grid),len(grid))

#treat command-line parameters 
mol       = read_from_file(sys.argv[1],ff=sys.argv[2])
coloumns  = int(sys.argv[3])
#example: 40,20,20
#these are only the counts in one direction
np_counts = np.array(map(int,sys.argv[4].split(",")))
#example: 0.35,0.5,0.5
np_del    = np.array(map(float,sys.argv[5].split(",")))
np_org    = np.array([0,0,0])

aligned_file = sys.argv[6]
dx_file      = sys.argv[7]
maxval       =float(sys.argv[8])

if len(sys.argv)>=10:
    correct = bool(sys.argv[9])
else:
    report=False
if len(sys.argv)>=11:
    report = bool(sys.argv[10])
else:
    report=False

#create helper vectors
np_start = np_org - np_counts*np_del
np_end   = np_org + np_counts*np_del

#this now contains the total number of counts per dimension
np_counts *= 2
np_counts += 1

#align the molecule's longest axis with the x-axis and the second longest axis with the y-direction
#and center the molecule to the origin
mol.align(list(np_org),[1,0,0],[0,1,0])
mol.write(aligned_file)
#append the same molecule again without translating or rotating it
mol.append(mol)
#since python needs quite some time to access an objects member, saving a member saves time
obmol = mol.mol
obff  = mol.ff
#create grid for translations of the molecule
space    = [np.linspace(start,end,num=counts,dtype=float) for start,end,counts in zip(np_start,np_end,np_counts)]
x,y,z    = np.array(np.meshgrid(*space,indexing="ij"))
x.shape  = (-1,1)
y.shape  = (-1,1)
z.shape  = (-1,1)
np_grid  = np.concatenate((x,y,z),axis=1)
grid     = double_dist(np_grid)

reset_vec = grid[len(grid)-1]
grid = [e for e in grid[:len(grid)-1]]

energies = list(gen_trans_en(obmol,obff,grid,reset_vec,maxval,report=report))

if correct:
    try:
        actualmax = max((e for e in energies if not e==maxval))
    except ValueError:
        actualmax = maxval
    energies = [actualmax if e==maxval else e for e in energies]

delx = [np_del[0],0.0,0.0]
dely = [0.0,np_del[1],0.0]
delz = [0.0,0.0,np_del[2]]

print_dx_file(dx_file,np_counts,np_start,delx,dely,delz,energies,comment="angles=(0.0 , 0.0 , 0.0)")
