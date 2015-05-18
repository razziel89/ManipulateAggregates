"""
A useful collection of functions to read in different data files
"""

import sys
import numpy as np

def read_xyz(file):
    """
    Read in an xyz-file. The first argument to be returned will be a list
    of element names and the second will be a list of the Cartesian
    coordinates.

    file: the path to the file from which data is to be read in
    """
    f=open(file,'r')

    #read the lines in the given file into the variable lines
    #and remove the trailing newline characters by using .rstrip()
    lines = np.array([line.rstrip() for line in f])
    #close the file descriptor
    f.close()
    
    #try to get the number of atoms in the molecule
    #if this does not succeed, the file is probably not a valid
    #xyz-file
    try:
        nr_atoms=int(lines[0])
    except ValueError:
        raise ValueError("This is probably not a valid xyz-file since the first line does not contain an integer.")
    
    #the first two lines of an xyz file are not necessary, hence, they are removed
    lines=np.array(lines[2:])
    
    #for every line take the last three coloumns as co-ordinates for the atom
    #line.split() yields a list of the coloumns in line (whitespace separation)
    #of those, take only the last ones line.split()[1:]
    #map(float,L) gives a list whose elements are those of the list L but converted to floats
    coordinates=np.array([map(float,line.split()[1:]) for line in lines])
    #for every line take the element in the first coloumn as the name of the atom's element
    names=[line.split()[0] for line in lines]
    #if there are more entries in the file than specified in the header, ignore the additional entries
    if len(coordinates)>nr_atoms:
        coordinates=coordinates[:nr_atoms]
        names=names[:nr_atoms]
    return names,coordinates

def read_afm(file,zscale=1):
    """
    read in a file as obtained from afm measurements
    the z-coordinate can be scaled by zscale
    x- and y-coordinates will be multiples of 1 and all the data
    will be centered around 0,0
    """
    maximum=sys.float_info.min
    minimum=sys.float_info.max
    f=open(file,"r")
    points=[]
    (x,y)=(0.0,0.0)
    for line in f:
        if not line.startswith("#"):
            x=0.0
            for entry in line.rstrip().split():
                value=float(entry)
                if value>maximum:
                    maximum=value
                if value<minimum:
                    minimum=value
                p=[x,y,value]
                points.append(p)
                x=x+1.0
            y=y+1.0
    x-=1
    y-=1
    for p in points:
        p[2]=((p[2]-minimum)/(maximum-minimum))*zscale
        p[1]=p[1]-y/2.0
        p[0]=p[0]-x/2.0
        
    return points

