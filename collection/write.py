"""
A handy collection of functions to write different filetypes.
"""

import re

class CommentError(Exception):
    pass

def print_xyz(filename,names,coordinates,width="10.6",comment=None):
    """
    Write data to an xyz-file. A comment can be added. If no comment is given,
    the filename will be taken as the comment.

    filename: The name of the file (will be overwritten if it exists already)
              or an already existing file descriptor. This way, you can print to special
              files like stdout or stderr.  Use sys.stdout or sys.stderr for this
              purpose.
    names: A list of strings containing the names of the atoms.
    coordinates: A list of 3-element lists containing the cartesian coordinates.
    comment: The content of the comment line as one string.
             Do not use newline characters.
    """
    if isinstance(filename, basestring):
        f=open(filename,'w')
        name=filename
    elif isinstance(filename, file):
        f=filename
        name=f.name
    else:
        raise TypeError("Specified file is neither a file descriptor nor a filename.")
    if comment==None:
        comment=name
    else:
        if re.search("\n",comment)!=None:
            raise CommentError("Specified comment contains a newline, which is not supported.")
    f.write(str(len(names))+"\n"+comment+"\n")
    for i in range(0,len(names)):
        tempstring=("%s    %"+width+"f    %"+width+"f    %"+width+"f\n")%(names[i],coordinates[i][0],coordinates[i][1],coordinates[i][2])
        f.write(tempstring)
    if isinstance(filename, basestring):
        f.close()

def _gen_cols(data,cols):
    i=0
    l=len(data)/cols
    while i<l:
        yield [data[cols*i+j] for j in xrange(3)]
        i+=1
    yield data[cols*i:]


def print_dx_file(filename,counts_xyz,org_xyz,delta_x,delta_y,delta_z,data,coloumns=3,comment=None):

    if isinstance(filename, basestring):
        f=open(filename,'w')
        name=filename
    elif isinstance(filename, file):
        f=filename
        name=f.name
    else:
        raise TypeError("Specified file is neither a file descriptor nor a filename.")
    if comment is None:
        comment="#"+name
    else:
        if re.search("\n",comment)!=None:
            raise CommentError("Specified comment contains a newline, which is not supported.")
        if not comment.startswith("#"):
            comment = "#"+comment
        if not comment.endswith("\n"):
            comment = comment+"\n"
        f.write(comment)
    #write header
    f.write("object 1 class gridpositions counts %4i %4i %4i\n"%tuple(counts_xyz))
    f.write("origin %7.6e %7.6e %7.6e\n"%tuple(org_xyz))
    f.write("delta %7.6e %7.6e %7.6e\n"%tuple(delta_x))
    f.write("delta %7.6e %7.6e %7.6e\n"%tuple(delta_y))
    f.write("delta %7.6e %7.6e %7.6e\n"%tuple(delta_z))
    f.write("object 2 class gridconnections counts %4i %4i %4i\n"%tuple(counts_xyz))
    prod=1
    for p in counts_xyz:
        prod *= p
    f.write("object 3 class array type double rank 0 items %12i data follows\n"%(prod))

    #write data
    for entry in _gen_cols(data,coloumns):
        tmp="%7.6e "*len(entry)+"\n"
        f.write(tmp%tuple(entry))
    
    if isinstance(filename, basestring):
        f.close()
