"""A handy collection of functions to write different filetypes.

@package ManipulateAggregates.collection.write
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

import re

class CommentError(Exception):
    """Raised if a comment is not valid."""
    pass

def print_xyz(filename,names,coordinates,width="10.6",comment=None):
    """Write data to an xyz-file.
    
    A comment can be added. If no comment is given, the filename will be taken
    as the comment.

    Args:
        filename: (string) the name of the file (will be overwritten if it
            exists already) or an already existing file descriptor. This way,
            you can print to special files like stdout or stderr. Use
            sys.stdout or sys.stderr for this purpose.
        names: (list of strings) a list of strings containing the names of the atoms.
        coordinates: (list of 3-element lists) contains the cartesian coordinates.

    Kwargs:
        width: (string) a format string that will be used to convert floats to
            strings. Defaults to "10.6".
        comment: (string) The content of the comment line as one string.
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
    if comment is None:
        comment=name
    else:
        if re.search(r"\n",comment)!=None:
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

def print_dx_file(filename,counts_xyz,org_xyz,delta_x,delta_y,delta_z,data,coloumns=3,comment=None, gzipped=False):
    """Print a dx file.

    Args:
        filename: (string) the file name
        counts_xyz: (tuple of 3 ints) how many points in each of the 3
            directions there are to the volumetric data
        org_xyz: (tuple of 3 floats) the origin of the volumetric data
        delta_x: (tuple of 3 floats) the Cartesian direction of the first
            voxel vector
        delta_y: (tuple of 3 floats) the Cartesian direction of the second
            voxel vector
        delta_z: (tuple of 3 floats) the Cartesian direction of the third
            voxel vector
        data: (list of floats) the volumetric data 

    Kwargs:
        coloumns: (int) in how many coloumns the volumetric data shall be
            written to the file
        comment: (string) a comment that is added at the top of the file
        gzipped: (bool) whether or not to write the file in gzipped format
    """
    if isinstance(filename, basestring):
        name=filename
        if gzipped:
            try:
                ##Try to import gzip. This wil throw an ImportError if that cannot be done.
                #import gzip
                #f=gzip.open(filename,"wb")
                from subprocess import Popen, PIPE
                process = Popen(['gzip', '--fast', '-c', '-'], stdin=PIPE, stdout=open(filename,'wb'), bufsize=4096)
                f = process.stdin
            except ImportError:
                print >>sys.stderr,"WARNING: cannot import gzip module, will treat %s as a non-gzipped one."%(filename)
                gzipped=False
                f=open(filename,"wb")
            except OSError:
                print >>sys.stderr,"WARNING: cannot import gzip module, will treat %s as a non-gzipped one."%(filename)
                gzipped=False
                f=open(filename,"wb")
        else:
            f=open(filename,"wb")
    elif isinstance(filename, file):
        f=filename
        name=f.name
        gzipped=False
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
        if len(entry)>0:
            f.write(tmp%tuple(entry))
    
    #write footer
    f.write('attribute "dep" string "positions"\n')
    f.write('object "regular positions regular connections" class field\n')
    f.write('component "positions" value 1\n')
    f.write('component "connections" value 2\n')
    f.write('component "data" value 3\n')

    if isinstance(filename, basestring):
        f.close()
        if gzipped:
            process.wait()

from openbabel import etab #used to transform element names and numbers
def _line_from_element_name(element,count,x,y,z):
    return "%s   %d   %d     %8.5f   %8.5f   %8.5f\n"%(element,count,etab.GetAtomicNum(element),x,y,z)

def _line_from_element_number(element,count,x,y,z):
    return "%s   %d   %d     %8.5f   %8.5f   %8.5f\n"%(etab.GetSymbol(element),count,element,x,y,z)

def _orbital_section(orb,count):
    result="%5d %d\n"%(count,0)
    for shell in orb[1]:
        orbtype,prefactor,nr_primitives,primitive_data = shell
        result+=" %s    %d %.2f\n"%(orbtype,nr_primitives,prefactor)
        for exponent,coefficient in primitive_data:
            result+="            %.6f    %.6f\n"%(exponent,coefficient)
    result+="\n"
    return result

def print_molden(filename,positions=None,pos_unit_string="Angs",element_names=True,GTO=None,MO=None):
    """Print a molden file.

    Args:
        filename: (string) the name of the file in which to save the data.

    Kwargs:
        positions: (list of [int,[float,float,float]] or list of [str,[float,float,float]])
            Contains information about atomic positions. The first entry
            defines the atom via a string or an int. @a element_names determines
            which one has been provided.
        pos_unit_string: (string) the string that will be put at the position
            where a programme expects the unit declaration. Some programmes
            seem to expect Bohr, some others AU or (AU).
        element_names: (bool) if True, positions has to be a list of
            [int,[float,float,float]]. Otherwise it has to be a list of
            [char,[float,float,float]]. Input for this can be generated by the
            function @a ManipulateAggregates.collection.read.molden_positions
        GTO: (list of [int, [[ str,float,int,[[float,float],[float,float],...], ...] ])
            Contains information about the Gaussian type orbital data. The
            first int specifies the number of shells in this orbital. P-type
            counts three times, F-type counts six times, etc. The first str
            declares the type of shell. The first float declares a general
            scaling factor for this orbital. The second int declares the number
            of prmitives in this shell. The [float,float] constructs define a
            primitive each as [exponent,prefactor]. Three dots indicate that
            the previous item can be repeated. Input for this can be generated
            by @a ManipulateAggregates.collection.read.molden_GTO
        MO: (list of [float,str,[float,...]]) Contains information about the
            molecular orbital coefficients. The first float declares the energy
            of the MO, the first str declares the spin ('alpha' or 'beta').
            The following list [float,...] contains one coefficient per shell
            (i.e. one for each S-type shell, 3 for each P-type shell, 6 for
            each F-type shell, etc.). Input for this can be generated using the
            function @a ManipulateAggregates.collection.read.molden_MO
    """
    if isinstance(filename, basestring):
        f=open(filename,'w')
        name=filename
    elif isinstance(filename, file):
        f=filename
        name=f.name
    else:
        raise TypeError("Specified file is neither a file descriptor nor a filename.")

    #write header
    f.write("[Molden Format]\n[Title]\nWritten by FireDeamon\n")

    #write atom positions if data has been provided
    #data has to be in Angstroms
    if positions is not None:
        f.write("[Atoms] %s\n"%(pos_unit_string))
        if element_names:
            linefunc=_line_from_element_name
        else:
            linefunc=_line_from_element_number
        count=1
        for element,(x,y,z) in positions:
            f.write(linefunc(element,count,x,y,z))
            count+=1

    #write GTO section if it has been provided
    if GTO is not None:
        f.write("[GTO]\n")
        count=1
        for orb in GTO:
            f.write(_orbital_section(orb,count))
            count+=1

    #write the MO section if it has been provided
    if MO is not None:
        f.write("[MO]\n")
        for orb in MO:
            energy,spin,occupation,coefficients = orb
            f.write(" Ene= %10.4f\n Spin= %s\n Occup= %.1f\n"%(energy,spin,occupation))
            count=1;
            for c in coefficients:
                f.write(" %5d %10.5f\n"%(count,c))
                count+=1

    if isinstance(filename, basestring):
        f.close()
