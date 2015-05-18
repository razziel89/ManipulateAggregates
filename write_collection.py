"""
A handy collection of functions to write different filetypes.
"""

def print_xyz(filename,names,coordinates,width="10.6"):
    """
    Write data to an xyz-file. The filename will be the comment line.
    filename: the name of the file (will be overwritten if it exists already)
    names: a list of strings containing the names of the atoms
    coordinates: a list of 3-element lists containing the cartesian coordinates
    """
    f=open(filename,'w')
    f.write(str(len(names))+"\n"+filename+"\n")
    for i in range(0,len(names)):
        tempstring=("%s    %"+width+"f    %"+width+"f    %"+width+"f\n")%(names[i],coordinates[i][0],coordinates[i][1],coordinates[i][2])
        f.write(tempstring)
    f.close()

