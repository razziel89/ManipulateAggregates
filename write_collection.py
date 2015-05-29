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
