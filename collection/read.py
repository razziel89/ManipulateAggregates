"""
A useful collection of functions to read in different data files
"""
import re
import sys
import numpy as np

class ReadCollectionError(Exception):
    pass

class WrongFormatError(ReadCollectionError):
    pass

class MissingSectionError(ReadCollectionError):
    pass

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
    #also ignore the last lines if there are more than the first line specifies
    lines=np.array(lines[2:nr_atoms+2])
    
    #for every line take the last three coloumns as co-ordinates for the atom
    #line.split() yields a list of the coloumns in line (whitespace separation)
    #of those, take only the last ones line.split()[1:]
    #map(float,L) gives a list whose elements are those of the list L but converted to floats
    coordinates=np.array([map(float,line.split()[1:]) for line in lines])
    #for every line take the element in the first coloumn as the name of the atom's element
    names=[line.split()[0] for line in lines]
    #if there are more entries in the file than specified in the header, ignore the additional entries
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

def _molden_positions(f,convert,elementnames=True,regex="^\[.+\]"):
    coords=[]
    line=f.next().rstrip()
    try:
        while not re.match(regex,line):
            l=line.split()
            if elementnames:
                at_name=l[0]
            else:
                at_name=int(l[2])
            at_pos=map(float,l[3:6])
            coords.append([at_name,[a*convert for a in at_pos]])
            line=f.next().rstrip()
    except (IndexError,ValueError) as e:
        raise WrongFormatError("Error on the current line: "+line+" ERROR is as follows: ",e)
    except StopIteration:
        pass
    return line,coords

def _molden_GTO(f,GTO_coefficients=False,nr_primitives=True,regex="^\[.+\]",int_regex="^(-|)[0-9]+$",float_regex="^(-|)[0-9]+(\.[0-9]*){0,1}$",orbital_regex="^(s|p|d|f|g|S|P|D|F|G)$"):
    """
    Read GTO section of a molden file.

    f: a handle to the molden file
    GTO_coefficients: same as for read_molden
    nr_primitives: same as GTO_nr_primitives for read_molden
    regex: a regular expression declaring the start of any section
           The section is considered to end either and EOF or
           when this regex matches the current line.
    int_regex: a regex matching an integer number
    float_regex: same for floats
    orbital_regex: a regex matching all valid orbital names
    """
    dict_orb_numbers={"s":1,"p":3,"d":6,"f":10,"g":15}
    gto=[]
    at_gto=0
    if GTO_coefficients:
        at_gto_coeff=[]
        shell=[]
    at_nr=0
    count=0
    elements=0
    line=f.next().rstrip()
    try:
        while not re.match(regex,line):
            l=line.split()
            if len(l)==2 and re.match(int_regex,l[0]) and re.match(int_regex,l[1]):
                if GTO_coefficients or nr_primitives:
                    if not count==elements:
                        raise WrongFormatError("Atom "+str(at_nr)+" needs "+str(elements)+" elements for basis function but only "+str(count)+" are available.")
                if at_gto>0 and at_nr>0:
                    if GTO_coefficients:
                        at_gto_coeff.append([orbitaltype,stretch,elements,shell])
                        shell=[]
                        gto.append([at_gto,at_gto_coeff])
                        at_gto_coeff=[]
                    else:
                        gto.append(at_gto)
                at_gto=0
                at_nr+=1
            elif len(l)==3 and re.match(orbital_regex,l[0]) and re.match(int_regex,l[1]) and re.match(float_regex,l[2]):
                if GTO_coefficients:
                    if len(shell)>0:
                        at_gto_coeff.append([orbitaltype,stretch,elements,shell])
                    shell=[]
                orbitaltype=l[0]
                stretch=float(l[2])
                if GTO_coefficients or nr_primitives:
                    elements=int(l[1])
                    count=0
                if nr_primitives:
                    at_gto+=dict_orb_numbers[l[0].lower()]*elements
                else:
                    at_gto+=dict_orb_numbers[l[0].lower()]
            elif len(l)==2 and re.match(float_regex,l[0]) and re.match(float_regex,l[1]):
                if nr_primitives or GTO_coefficients:
                    count+=1
                if GTO_coefficients:
                    shell.append(map(float,[l[0],l[1]]))
            elif re.match("^\s*$",line):
                pass
            line=f.next().rstrip()
    except (IndexError,ValueError) as e:
        raise WrongFormatError("Error on the current line: "+line+" ERROR is as follows: ",e)
    except StopIteration:
        pass
    if at_nr>0:
        if GTO_coefficients:
            at_gto_coeff.append([orbitaltype,stretch,elements,shell])
            sheel=[]
            gto.append([at_gto,at_gto_coeff])
        else:
            gto.append(at_gto)
    return line,gto

def _molden_MO(f,MO_coefficients=False,regex="^\[.+\]",int_regex="^(-|)[0-9]+$",float_regex="^(-|)[0-9]+(\.[0-9]*){0,1}$"):
    """
    Read MO section of a molden file.

    f: a handle to the molden file
    regex: a regular expression declaring the start of any section
           The section is considered to end either and EOF or
           when this regex matches the current line.
    int_regex: a regex matching an integer number
    float_regex: same for floats
    """
    mo=[]
    orbital=[]
    energy=None
    #is_alpha==True means it is alpha, False means it's beta
    is_alpha=None
    occupation=None
    line=f.next().rstrip()
    try:
        while not re.match(regex,line):
            l=line.split()
            if len(l)==2 and re.match("^.+=$",l[0]):
                if MO_coefficients:
                    if len(orbital)>0:
                        if not occupation==None and not energy==None and not is_alpha==None:
                            mo.append([energy,"alpha" if is_alpha else "beta",occupation,orbital])
                            energy=None
                            is_alpha=None
                            occupation=None
                        else:
                            raise WrongFormatError("At least one orbital does not specify all of Ene=, Spin= and Occup=")
                elif not occupation==None and not energy==None and not is_alpha==None:
                    mo.append([energy,"alpha" if is_alpha else "beta",occupation])
                    energy=None
                    is_alpha=None
                    occupation=None
                if MO_coefficients:
                    orbital=[]
                    count=0
                if l[0]=="Ene=":
                    energy=float(l[1])
                elif l[0]=="Spin=":
                    if re.match("^[Aa][Ll][Pp][Hh][Aa]$",l[1]):
                        is_alpha=True
                    elif re.match("^[Bb][Ee][Tt][Aa]$",l[1]):
                        is_alpha=False
                    else:
                        raise WrongFormatError("Spin name has to be either alpha or beta but it is "+l[1])
                elif l[0]=="Occup=":
                    occupation=float(l[1])
            elif len(l)==2 and re.match(int_regex,l[0]) and re.match(float_regex,l[1]):
                if MO_coefficients:
                    newcount=int(l[0])
                    #if contributions have been left out from the molden file, fill them with zeroes
                    for i in xrange(count+1,newcount):
                        orbital.append(0.0)
                    count=newcount
                    orbital.append(float(l[1]))
            elif re.match("^\s*$",line):
                pass
            line=f.next().rstrip()
    except (IndexError,ValueError,WrongFormatError) as e:
        raise WrongFormatError("Error on the current line: "+line+" ERROR is as follows: ",e)
    except StopIteration:
        pass
    if MO_coefficients:
        if len(orbital)>0:
            if not occupation==None and not energy==None and not is_alpha==None:
                mo.append([energy,"alpha" if is_alpha else "beta",occupation,orbital])
            else:
                raise WrongFormatError("At least one orbital does not specify all of Ene=, Spin= and Occup=")
    elif not occupation==None and not energy==None and not is_alpha==None:
        mo.append([energy,"alpha" if is_alpha else "beta",occupation])

    if MO_coefficients:
        max_nr_primitives=max([len(o[-1]) for o in mo])
        #set all orbitals that have been filled to the same maximum range
        for i in xrange(len(mo)):
            mo[i][-1]+=[0.0]*(max_nr_primitives-len(mo[i][-1]))
    return line,mo

def read_molden(file,positions=True,elementnames=True,GTO=True,GTO_coefficients=False,GTO_nr_primitives=False,MO=True,MO_coefficients=False):
    """
    This function read in a Molden file accordig to some flags that are st.
    Beware: only Cartiesian coordinates are supported so far. Will return a
    dictionary with approrpiately named entries. Only Cartersian Gaussian type
    orbitals are supported as of now.

    file: the path to the file from which data is to be read in
    positions: whether or not atomic coordinates shall be read in. Will
               always be returned in Angstroms.
    elementnames: whether elements shall be identified by their names (True)
                  or by their element numbers (False)
    GTO: whether or not to read in the GTO section. True will result in
         a list assigning orbitals to a particular atom.
    GTO_coefficients: whether or not to read in all GTO-coefficients
                      If False, only the number of primitives or shells
                      are counted and returned (influenced by GTO_nr_primitives)
    GTO_nr_primitives: whether you want to count the number of primitives (True)
                       or shells (False)
    MO: whether or not to read in the MO-section. You will get 2 lists
        one for each spin containing energy and occupation
    MO_coefficients: whether or not all MO-coefficients shall be read in
                     If False, only auxilliary information line energies and
                     occupations and spins are read in.
    """
    #regular expression that matches the beginning of any secion
    regex="^\[.+\]"
    #regular expression that matches an integer
    int_regex="^(-|)[0-9]+$"
    #regular expression that matches a float
    float_regex="^(-|)[0-9]+(\.[0-9]*){0,1}$"
    #regular expression that matches any of the allowed orbital names
    orbital_regex="^(s|p|d|f|g|S|P|D|F|G)$"
    #open the file
    f=open(file,"r")
    #initialize dictionary that will hold the results
    result={}
    #check format of first line
    if not re.match("^\[Molden Format\]\s*$",f.next().rstrip()):
        raise WrongFormatError("The first line in a Molden file has to be '[Molden Format]'")
    nr_sections=[positions,GTO,MO].count(True)
    sec=nr_sections-1
    important_sections_regex=""
    if positions:
        important_sections_regex+="^\[Atoms\]\s+([Aa][Nn][Gg][Ss]|[Aa][Uu])\s*$"
        if sec>0:
            important_sections_regex+="|"
            sec-=1
    if GTO:
        important_sections_regex+="^\[GTO\]"
        if sec>0:
            important_sections_regex+="|"
            sec-=1
    if MO:
        important_sections_regex+="^\[MO\]"
        if sec>0:
            important_sections_regex+="|"
            sec-=1
    sec=0
    try:
        #skip to the next important section
        line=f.next().rstrip()
        while not re.match(important_sections_regex,line):
            line=f.next().rstrip()
        #read in only the requested number of sections
        while sec<nr_sections:
            #read in requested sections
            if re.match("^\[Atoms\]\s+([Aa][Nn][Gg][Ss]|[Aa][Uu])\s*$",line) and positions:
                #Atoms section
                sec+=1
                #determine whether the atomic coordinates are giveb in angstroms or bohrs
                if re.match("^\[Atoms\]\s+([Aa][Nn][Gg][Ss])\s*$",line):
                    convert=1
                else:
                    convert=0.52918
                #read in the section
                line,result["positions"] = _molden_positions(f,convert,elementnames=elementnames,regex=regex)
            elif re.match("^\[GTO\]",line) and GTO:
                #GTO section
                sec+=1
                #read in the section
                line,result["GTO"] = _molden_GTO(f,GTO_coefficients=GTO_coefficients,regex=regex,nr_primitives=GTO_nr_primitives,int_regex=int_regex,float_regex=float_regex,orbital_regex=orbital_regex)
            elif re.match("^\[MO\]",line) and MO:
                #MO section
                #GTO section
                sec+=1
                #read in the section
                line,result["MO"] = _molden_MO(f,MO_coefficients=MO_coefficients,regex=regex,int_regex=int_regex,float_regex=float_regex)

            #break loop prematurely if enough sections have been read in
            if re.match(regex,line):
                if sec>=nr_sections:
                    break
            #skip to next significant section if it has not yet been reached
            if not re.match(important_sections_regex,line):
                while not re.match(important_sections_regex,line):
                    line=f.next().rstrip()
    except StopIteration:
        if sec<nr_sections:
            #if end of file reached before the requested sections could be read in
            raise WrongFormatError("You requested "+str(nr_sections)+" sections but only "+str(sec)+" were found as end of file was reached.")
    #do some sanity checks on the results
    if positions and len(result["positions"])==0:
        raise MissingSectionError("Atomic coordinates requested but whole file read in without finding the secion.")
    if GTO and len(result["GTO"])==0:
        raise MissingSectionError("GTO section requested but whole file read in without finding the secion.")
    if MO and len(result["MO"])==0:
        raise MissingSectionError("MO section requested but whole file read in without finding the secion.")
    f.close()
    return result
