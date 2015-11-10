"""
A useful collection of functions to read in different data files
"""
import re
import sys
import numpy as np
import itertools
import ConfigParser

class ReadCollectionError(Exception):
    pass

class WrongFormatError(ReadCollectionError):
    pass

class MissingArgumentError(ReadCollectionError):
    pass

class MissingSectionError(ReadCollectionError):
    pass

class MoleculesNotMatchingError(ReadCollectionError):
    pass

class CannotMatchError(ReadCollectionError):
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
                    tempshell = list(map(float,[l[0],l[1]]))
                    shell.append(tempshell)
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

    All units will be converted to Bohr when reading in the file. This
    function assumes that all exponents in the GTO section are given in Bohr^-1

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
                #determine whether the atomic coordinates are given in angstroms or bohrs
                if re.match("^\[Atoms\]\s+([Aa][Nn][Gg][Ss])\s*$",line):
                    convert=1.0/0.52918
                else:
                    convert=1.0
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

def _renormalize_whole(vec,norm=1.0):
    return vec*(norm/np.linalg.norm(vec))

def _renormalize_individual(vec,norm=1.0):
    return vec*(norm/np.max(np.linalg.norm(vec,axis=1)))

def _renormalize_none(vec,norm=1.0):
    return vec*norm

def read_aims_frequencies(fname,mode=1,amplitude_factor=1,normalize="individual"):
    """
    Read in a frequencies file in AIMS format.

    fname: filename
    mode: mode to be read in (starting at 1)
    amplitude_factor: the new norm of the vector according to the value of normalize
    normalize: possible values: whole, individual, none
               whole:      Normalize the whole displacement vector
               individual: For the atom with the largest total displacement vector, normalize
                           this to amplitude_factor and adjust all the others accordingly.
               none:       Do not perform normalization, only scale by amplitude_factor
    """
    #read the lines in the given file into the variable lines
    #and remove the trailing newline characters by using .rstrip()
    f=open(fname)
    line1 = f.next().rstrip() 
    
    nr_atoms=int(line1.split()[0])
    
    nr_deg_of_freedom=3*nr_atoms
    nr_normalmodes=int(line1.split()[1])
    
    displacement=np.zeros(nr_deg_of_freedom)
    
    mode_count=0
    while mode_count <= mode:
    	line = f.next().rstrip()
    	freqs=map(float,line.split())
    	for value in freqs:
    		mode_count+=1
    		if mode_count == mode:
    			frequency=value
    while mode_count < nr_normalmodes:
    	mode_count += len(f.next().rstrip().split())
    
    coord=0
    while coord < mode*nr_deg_of_freedom:
    	line = f.next().rstrip()
    	disp=map(float,line.split())
    	for value in disp:
    		if coord>=(mode-1)*nr_deg_of_freedom and coord<mode*nr_deg_of_freedom:
    			displacement[coord%nr_deg_of_freedom]=value
    		coord+=1
    f.close()
    displacement.shape=(nr_atoms,3)
    if normalize=="individual":
        displacement=_renormalize_individual(displacement,amplitude_factor)
    elif normalize=="whole":
        displacement=_renormalize_whole(displacement,amplitude_factor)
    else:
        displacement=_renormalize_none(displacement,amplitude_factor)

    return frequency,displacement

def read_terachem_frequencies(fname,mode=1,amplitude_factor=1,normalize="individual"):
    """
    Read in a frequencies file in TeraChem format.

    fname: filename
    mode: mode to be read in (starting at 1)
    amplitude_factor: the new norm of the vector according to the value of normalize
    normalize: possible values: whole, individual, none
               whole:      Normalize the whole displacement vector
               individual: For the atom with the largest total displacement vector, normalize
                           this to amplitude_factor and adjust all the others accordingly.
               none:       Do not perform normalization, only scale by amplitude_factor
    """
    f=open(fname)
    #read the lines in the given file into the variable lines
    #and remove the trailing newline characters by using .rstrip()
    line1 = f.next().rstrip() 
    
    nr_atoms=int(line1.split()[2])
    
    nr_deg_of_freedom=3*nr_atoms

    #number of normalmodes is in second line
    line1 = f.next().rstrip() 
    nr_normalmodes=int(line1.split()[3])

    #the third line does not contain any useful information
    f.next()
    
    displacement=np.zeros(nr_deg_of_freedom)
    disp_count=0

    #mode counting starts at 0
    mode-=1
    maxmode=-1

    while maxmode<mode:
        modes=map(int,f.next().rstrip().split())
        freqs=map(float,f.next().rstrip().split())
        #the next line does not contain any useful information
        f.next()
        maxmode=max(modes)
        if maxmode<mode:
            for i in range(nr_atoms*3+1):
                f.next()
        else:
            index=modes.index(mode)
            frequency=freqs[index]
            for at in range(nr_atoms):
                #treat first line per atom
                entries=f.next().rstrip().split()
                atom=int(entries[0])
                disps=map(float,entries[1:])
                displacement[disp_count]=disps[index]
                disp_count+=1
                #treat second line per atom
                entries=f.next().rstrip().split()
                disps=map(float,entries)
                displacement[disp_count]=disps[index]
                disp_count+=1
                #treat third line per atom
                entries=f.next().rstrip().split()
                disps=map(float,entries)
                displacement[disp_count]=disps[index]
                disp_count+=1

    f.close()
    displacement.shape=(nr_atoms,3)
    if normalize=="individual":
        displacement=_renormalize_individual(displacement,amplitude_factor)
    elif normalize=="whole":
        displacement=_renormalize_whole(displacement,amplitude_factor)
    else:
        displacement=_renormalize_none(displacement,amplitude_factor)

    return frequency,displacement

def read_charges_simple(file,compare_elements=False,molecule=None):
    """
    Read in an xyz-file where each line of Cartesian coordinates is followed by
    a charge. The first argument to be returned will be the position of the
    partial charges and the second will be a list of the charges.

    file: the path to the file from which data is to be read in
    compare_elements: if this is True, molecule must be a molecule object.
                      A sanity check will be performed where the element
                      names from the molecule object are compared to those
                      from the given file. Furthermore, the coordinates
                      are taken not from the file but from the molecule
                      object.
    molecule: the molecule object to compare against
    """
    f=open(file,'r')

    #read the lines in the given file into the variable lines
    #and remove the trailing newline characters by using .rstrip()
    lines = np.array([line.rstrip().split() for line in f])
    #close the file descriptor
    f.close()
    
    #try to get the number of atoms in the molecule
    #if this does not succeed, the file is probably not a valid
    #xyz-file
    try:
        nr_atoms=int(lines[0][0])
    except ValueError:
        raise ValueError("This is probably not a valid xyz-file since the first line does not contain an integer.")
    
    #the first two lines of an xyz file are not necessary, hence, they are removed
    #also ignore the last lines if there are more than the first line specifies
    lines=np.array(lines[2:nr_atoms+2])
    
    if compare_elements:
        elements_molecule=[e for e in sorted(molecule.get_names())]
        elements_file=[line[0] for line in sorted(lines)]
        if elements_molecule==elements_file:
            charges=np.array([float(line[4]) for line in lines])
            coordinates=np.array(molecule.get_coordinates())
        else:
            raise MoleculesNotMatchingError("Molecule read from the charge file and the given molecule object do not contain the same elements.")
    else:
        try:
            coordinates=np.array([map(float,line[1:4]) for line in lines])
            charges=np.array([float(line[4]) for line in lines])
        except IndexError:
            raise IndexError("Not enough coloumns! There need to be 4 in every line but the first to. Element name, x,y,z coordinates, charge.")
        except ValueError:
            raise ValueError("At least one value on one of the lines is no valid float.")
    return coordinates, charges

def _list_equiv(l1,l2):
    """
    Check if two lists are equivalent by comparing them element wise.
    If one list is exhausted and they were equivalent so far, consider
    them equivalent.
    """
    for e1,e2 in zip(l1,l2):
        if not e1==e2:
            return False
    return True

def read_charges_dx(file,add_nuclear_charges=False,molecule=None,unit_conversion=1.0,invert_charge_data=True,rescale_charges=True,total_charge=0,nr_return=None,density=False,header_dict=None):
    """
    Read in a DX file. Will return Cartesian coordinates of data
    and data.

    file:               the name of the dx file
    add_nuclear_charges:if True, nuclear charges and coordinates will be the first 
                        entries in the file, needs molecule to be declared
    molecule:           object of class molecule. get_coordinates() and get_charges()
                        methods are used to get the coordinates and values of nuclear
                        point charges.
    unit_conversion:    give a value by which to scale all coordinates (might displace center)
    invert_charge_data: if False, do not invert the volumetric charge data. Nuclear charges
                        are always positive and volumetric data is taken inverted by default.
    rescale_charges:    the sum of atomic charges and the sum of nuclear charges have to
                        match. If this is True and the charges don't match, rescale all 
                        volumetric data linearly so that they do. 
                        Only makes sense if add_nuclear_charges == True
    total_charge:       the total charge of the molecule to properly rescale the electronic
                        charges
    nr_return:          if a variable of type list is given, append to it the number of
                        atoms and the number of volumetric entries
    density:            if True, return the density at the center of the voxel instead of
                        the product of the density and the voxel's volume
    header_dict:        if given a dictionary, add to it the entries counts_xyz,org_xyz,delta_x,delta_y
                        and delta_z to allow easy re-printing of the dx-file
    """
    dxdata      = read_dx(file,unit_conversion=unit_conversion,invert_charge_data=invert_charge_data,density=density,header_dict=header_dict,grid=True,data=True)
    charges     = np.array(dxdata["data"])
    coordinates = np.array(dxdata["grid"])

    if add_nuclear_charges:
        if not molecule==None:
            try:
                #since no data about the nuclei is saved in the dx file,
                #get it from the molecule object
                nuc_charges=np.array(molecule.get_charges())
                nuc_coordinates=np.array(molecule.get_coordinates())

            except AttributeError as e:
                raise AttributeError("Given molecule object does not define methods get_charges or get_coordinates.",e)
        else:
            raise MissingArgumentError("Cannot add nuclear charges without molecule argument.")

        if rescale_charges:
            #This equatiom makes it so that the sum over the electronic charges (which is negative)
            #plus the sum of the nuclear charges (which is positive) equals the total charge.
            charges *= (total_charge - np.sum(nuc_charges)) / np.sum(charges)

        np.insert(coordinates,0,nuc_coordinates,axis=0)
        np.insert(charges,0,nuc_charges,axis=0)

    return coordinates,charges

def read_dx(file,unit_conversion=1.0,invert_charge_data=False,density=True,header_dict=None,grid=True,data=True,silent=False):
    """
    Read in a DX file. Will return a dictionary with the entries grid (if grid==True) and data (if data == True).

    unit_conversion:    give a value by which to scale all coordinates (might displace center)
    invert_charge_data: if False, do not invert the volumetric charge data.
    density:            if True, return the density at the center of the voxel instead of
                        the product of the density and the voxel's volume
    header_dict: empty dictionary
        If given a dictionary, add to it the entries counts_xyz,org_xyz,delta_x,delta_y
        and delta_z to allow easy re-printing of the dx-file.
    grid: boolean
        Whether or not to return the grid.
    data: boolean
        Whether or not to return the data.
    silent: boolean
        Whether or not to utter warnings such as about a missing footer.
    """
    result = {}
    f=open(file)
    if header_dict is not None:
        import copy
    #ignore comment lines at the beginning
    l=f.next()
    while l.startswith("#"):
        l=f.next()
    #according to the VMD mailing list, the ordering is: z fast, y medium, and x slow
    #this translates to: z inner, y middle, and x outer 
    #Source: http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/21526.html
    match_dict={'outer':0,'inner':2,'middle':1}
    xyz_dict={'x':0,'y':1,'z':2}
    #preallocate stuff
    axes=[None]*3
    nrs=[None]*3
    #Comment on units: the grid (if requested) will have the units used in the file itself multiplied
    #                  by unit_conversion
    #read in the header
    #line stating the number of gridpositions
    line=l.rstrip().split()
    if _list_equiv(line,["object","1","class","gridpositions","counts"]):
        try:
            nrs=map(int,line[5:])
        except ValueError as e:
            raise ValueError("First non-comment line in DX file does not end on three integers.",e)
    else:
        raise WrongFormatError("First non-comment line in DX file must be 'object 1 class gridpositions counts nx ny nz'")
    if header_dict is not None:
        header_dict["counts_xyz"]=copy.copy(list(nrs))
    #line with origin
    line=f.next().rstrip().split()
    if line[0]=="origin":
        try:
            origin=np.array(map(float,line[1:]))*unit_conversion
        except ValueError as e:
            raise ValueError("Second non-comment line in DX file does not end on three floats.",e)
    else:
        raise WrongFormatError("Second non-comment line in DX file must be 'origin ox oy oz'")
    if header_dict is not None:
        header_dict["org_xyz"]=copy.copy(list(origin))

    #the three lines with the voxel sizes
    for c in range(3):
        line=f.next().rstrip().split()
        if line[0]=="delta":
            try:
                axes[c]=map(float,line[1:4])
            except ValueError as e:
                raise ValueError("One of third to fifth non-comment lines in DX file does not end on three floats.",e)
        else:
            raise WrongFormatError("Third to fifth non-comment lines must be 'delta dx dy dz'")
    if header_dict is not None:
        header_dict["delta_x"]=copy.copy(list(axes[0]))
        header_dict["delta_y"]=copy.copy(list(axes[1]))
        header_dict["delta_z"]=copy.copy(list(axes[2]))
    #next line, which somewhat of a duplicate, the integers are ignored here
    line=f.next().rstrip().split()
    if _list_equiv(line,["object","2","class","gridconnections","counts"]):
        try:
            map(int,line[5:])
        except ValueError as e:
            raise ValueError("Sixth non-comment line in DX file does not end on three integers.",e)
    else:
        raise WrongFormatError("Sixth non-comment line in DX file must be 'object 2 class gridconnections counts nx ny nz'")
    #next line contains some test data to check whether format is correct
    line=f.next().rstrip().split()
    if _list_equiv(line,["object","3","class","array","type","double","rank","0","items"]) and _list_equiv(line[10:],["data","follows"]):
        try:
            if not nrs[0]*nrs[1]*nrs[2] == int(line[9]):
                raise WrongFormatError("Seventh non-comment line in DX file does not contain the correct number of volumetric data elements.")
        except ValueError as e:
            raise ValueError("Seventh non-comment line in DX file does not contain the total number of volumetric data elements",e)
    else:
        raise WrongFormatError("Seventh non-comment line in DX file must be 'object 3 class array type double rank 0 items nx*ny*nz data follows'")
    #convert to numpy arrays
    axes=np.array(axes)
    nrs=np.array(nrs)
    #this is the volume of one voxel
    #which will be used to convert charge density to charge
    #so it will seem as if there were a point charge at the center
    #of the voxel containing the whole charge inside that voxel
    if density:
        volume=1.0
    else:
        volume=np.linalg.det(unit_conversion*axes.transpose())
    #read in volumetric data if requested
    if data:
        charges=np.zeros((nrs[0]*nrs[1]*nrs[2]),dtype=float)
        count=0
        for l in f:
            if l.startswith(("0","1","2","3","4","5","6","7","8","9","+","-")):
                for e in map(float,l.rstrip().split()):
                    charges[count]=e*volume
                    count+=1
            else:
                break
        if invert_charge_data:
            charges*=-1
        result["data"] = charges
    else:
        for l in f:
            if not l.startswith(("0","1","2","3","4","5","6","7","8","9","+","-")):
                break
    try:
        #read in footer which is assumed to be of a certain format
        #this might break if a programme changes data assignments
        #in l is the first line of the footer
        while l.startswith(("0","1","2","3","4","5","6","7","8","9","+","-")):
            l=f.next()
        l=l.rstrip().split()
        if not _list_equiv(l,["attribute",'"dep"',"string",'"positions"']):
            raise WrongFormatError('First line of footer must be attribute \'"dep" string "positions"\'')
        l=f.next().rstrip().split()
        if not _list_equiv(l,["object",'"regular',"positions","regular",'connections"',"class","field"]):
            raise WrongFormatError('Second line of footer must be attribute \'object "regular positions regular connections" class field\'')
        l=f.next().rstrip().split()
        if not _list_equiv(l,["component",'"positions"',"value","1"]):
            raise WrongFormatError('Third line of footer must be attribute \'component "positions" value 1\'')
        l=f.next().rstrip().split()
        if not _list_equiv(l,["component",'"connections"',"value","2"]):
            raise WrongFormatError('Third line of footer must be attribute \'component "connections" value 2\'')
        l=f.next().rstrip().split()
        if not _list_equiv(l,["component",'"data"',"value","3"]):
            raise WrongFormatError('Third line of footer must be attribute \'component "data" value 3\'')
    except StopIteration:
        if not silent:
            print >>sys.stderr,"WARNING: no footer present in dx-file"
    #preallocate
    if grid:
        coordinates        = np.zeros((nrs[0]*nrs[1]*nrs[2],3),dtype=float)
        axes_rearranged    = np.zeros(axes.shape,dtype=float)
        axes_rearranged[0] = axes[match_dict['outer']]
        axes_rearranged[1] = axes[match_dict['middle']]
        axes_rearranged[2] = axes[match_dict['inner']]
        nrs_rearranged     = np.zeros(nrs.shape,dtype=int)
        nrs_rearranged[0]  = nrs[match_dict['outer']]
        nrs_rearranged[1]  = nrs[match_dict['middle']]
        nrs_rearranged[2]  = nrs[match_dict['inner']]
        #this is a fast variant to do the rearranging of the data. see read_charges_cube for the other variants
        result["grid"]     = [
                                (origin+np.dot(multi_indices,axes_rearranged))*unit_conversion 
                             for multi_indices in np.indices(nrs_rearranged).reshape(3,-1).T
                             ]
    f.close()
    return result

def read_charges_cube(file,match_word_order=False,match_axis_order=True,add_nuclear_charges=False,force_angstroms=False,invert_charge_data=False,rescale_charges=True,total_charge=0,nr_return=None,density=False):
    """
    Read in a Gaussian-Cube file. Will return Cartesian coordinates of charges
    and charges.

    file:               the name of the cube file
    match_word_order:   if True, try to find out the order of inner, outer and middle to
                        guess the correct volumetric data.
                        The words outer, inner and middle have to be pressent in
                        a certain order. Example:
                        OUTER X, INNER Y, MIDDLE Z will result in x, y, z being the
                        outer, inner and middle loops, respectively.
                        Per default, outer, middle and inner loop are x,y and z, respectively.
    match_axis_order:   if True, try to find out the order of X, Y and Z coordinates of
                        the volumetric data.
                        The letters X, Y and Z have to be present in a certain order.
                        Default: see match_word_order
    add_nuclear_charges:if True, nuclear charges and coordinates will be the first 
                        entries in the file
    force_angstroms:    if True, enforce everything to be considered to be in Angstroms
    invert_charge_data: if True, invert the volumetric charge data. Nuclear charges
                        are always positive and volumetric data is taken as is.
    rescale_charges:    the sum of atomic charges and the sum of nuclear charges have to
                        match. If this is True and the charges don't match, rescale all 
                        volumetric data linearly so that they do. 
                        Only makes sense if add_nuclear_charges == True
    total_charge:       the total charge of the molecule to properly rescale the electronic
                        charges
    nr_return:          if a variable of type list is given, append to it the number of
                        atoms and the number of volumetric entries
    density:            if True, return the density at the center of the voxel instead of
                        the product of the density and the voxel's volume
    """
    f=open(file)
    f.next()
    line=f.next().rstrip()
    #per default, outer, middle and inner loop are x,y and z, respectively
    matching_wordperm=['outer','middle','inner']
    matching_xyzperm=['x','y','z']
    match_dict={'outer':0,'inner':0,'middle':0}
    xyz_dict={'x':0,'y':1,'z':2}
    #try to find out the order of the words outer, inner and middle as well as the order of x, y and z
    if match_word_order:
        if None in (re.search(regex,line,re.IGNORECASE) for regex in ['\\bouter\\b','\\binner\\b','\\bmiddle\\b']):
            raise CannotMatchError("You requested to match the word order against the second line but I cannot find the words outer, inner, middle there.")
        else:
            matching_wordperm=[]
            for perm in itertools.permutations(['outer','inner','middle']):
                if re.match('.*\\b'+'\\b.*\\b'.join(perm)+'\\b.*',line,re.IGNORECASE):
                    matching_wordperm=perm
                    break
    if match_axis_order:
        if None in (re.search(regex,line,re.IGNORECASE) for regex in ['\\bx\\b','\\by\\b','\\bz\\b']):
            raise CannotMatchError("You requested to match the order against the second line but I cannot find the words x, y, z there.")
        else:
            matching_xyzperm=[]
            for perm in itertools.permutations(['x','y','z']):
                if re.match('.*\\b'+'\\b.*\\b'.join(perm)+'\\b.*',line,re.IGNORECASE):
                    matching_xyzperm=perm
                    break
    for word,letter in zip(matching_wordperm,matching_xyzperm):
        match_dict[word]=xyz_dict[letter]
    #preallocate stuff
    axes=[None]*3
    nrs=[None]*3
    #prepare array for unit conversion
    #default in file is atomic units but everything will be transformed to Angstroms
    unit_conversion=np.array([0.5291772488]*3,dtype=float)
    #read in the header
    #line with origin and number of atoms
    line=f.next().rstrip().split()
    nr_atoms=int(line[0])
    origin=np.array(map(float,line[1:4]))
    #the three lines with the voxel sizes and numbers of entries
    for c in range(3):
        line=f.next().rstrip().split()
        nrs[c]=abs(int(line[0]))
        if int(line[0])<0 or force_angstroms:
            unit_conversion[c]=1.0
        axes[c]=map(float,line[1:4])
    #the units of the origin have to be converted as well! doh...
    oring=origin*unit_conversion
    #convert to numpy arrays
    axes=np.array(axes)
    nrs=np.array(nrs)
    #this is the volume of one voxel
    #which will be used to convert charge density to charge
    #so it will seem as if there were a point charge at the center
    #of the voxel containing the whole charge inside that voxel
    if density:
        volume=1.0
    else:
        volume=np.linalg.det(unit_conversion*axes.transpose())
    #read in volumetric data
    if add_nuclear_charges:
        sum_nuclear_charges=total_charge
        charges=np.zeros((nrs[0]*nrs[1]*nrs[2]+nr_atoms),dtype=float)
        coordinates=np.zeros((nrs[0]*nrs[1]*nrs[2]+nr_atoms,3),dtype=float)
        for count in xrange(nr_atoms):
            line=f.next().rstrip().split()
            #nuclear charges have to have the opposite sign as electronic charges
            charges[count]=int(line[0])
            sum_nuclear_charges+=charges[count]
            #the first coloumn contains the atomic charge and the second is undefined
            #so the last 3 contain the information I need`
            coordinates[count]=map(float,line[2:5])
        count=nr_atoms
    else:
        #skip lines of atomic positions
        charges=np.zeros((nrs[0]*nrs[1]*nrs[2]),dtype=float)
        coordinates=np.zeros((nrs[0]*nrs[1]*nrs[2],3),dtype=float)
        count=0
        for i in xrange(nr_atoms):
            f.next()
    sum_electronic_charges=0
    for l in f:
        for e in map(float,l.rstrip().split()):
            charges[count]=e*volume
            sum_electronic_charges+=charges[count]
            count+=1
    if add_nuclear_charges and ( rescale_charges or invert_charge_data):
        is_nucleus=np.zeros(charges.shape,dtype=bool)
        is_nucleus[:nr_atoms]=np.ones((nr_atoms),dtype=bool)
    if add_nuclear_charges and rescale_charges:
        electronic_charge_rescale_factor=sum_nuclear_charges/sum_electronic_charges
        charges=charges*is_nucleus+charges*np.logical_not(is_nucleus)*electronic_charge_rescale_factor
    if invert_charge_data:
        if add_nuclear_charges:
            charges=charges*is_nucleus+charges*np.logical_not(is_nucleus)*(-1)
        else:
            charges*=-1
    axes_rearranged=np.zeros(axes.shape,dtype=float)
    axes_rearranged[0]=axes[match_dict['outer']]
    axes_rearranged[1]=axes[match_dict['middle']]
    axes_rearranged[2]=axes[match_dict['inner']]
    nrs_rearranged=np.zeros(nrs.shape,dtype=int)
    nrs_rearranged[0]=nrs[match_dict['outer']]
    nrs_rearranged[1]=nrs[match_dict['middle']]
    nrs_rearranged[2]=nrs[match_dict['inner']]
    unit_conversion_rearranged=np.zeros(unit_conversion.shape,dtype=float)
    unit_conversion_rearranged[0]=unit_conversion[match_dict['outer']]
    unit_conversion_rearranged[1]=unit_conversion[match_dict['middle']]
    unit_conversion_rearranged[2]=unit_conversion[match_dict['inner']]
    if add_nuclear_charges:
        count=nr_atoms
    else:
        count=0
    #o, m, i stand for outer, inner and middle, respectively
    #All the 4 variants do the exact same thing, but the last is approx. twice as fast as the first
    #Variant 1
    #for multi_indices in (np.array((o,m,i)) for o in xrange(nrs_rearranged[0]) for m in xrange(nrs_rearranged[1]) for i in xrange(nrs_rearranged[2])):
    #    position=(origin+np.dot(multi_indices,axes_rearranged))*unit_conversion_rearranged
    #    coordinates[count]=position
    #    count+=1
    #Variant 2
    #for o in xrange(nrs_rearranged[0]):
    #    for m in xrange(nrs_rearranged[1]):
    #        for i in xrange(nrs_rearranged[2]):
    #            position=(origin+np.dot(np.array((o,m,i)),axes_rearranged))*unit_conversion_rearranged
    #            coordinates[count]=position
    #            count+=1
    #Variant 3
    #for multi_indices in np.indices(nrs_rearranged).reshape(3,-1).T:
    #    position=(origin+np.dot(multi_indices,axes_rearranged))*unit_conversion_rearranged
    #    coordinates[count]=position
    #    count+=1
    #variant 4
    coordinates[count:]=[(origin+np.dot(multi_indices,axes_rearranged))*unit_conversion_rearranged for multi_indices in np.indices(nrs_rearranged).reshape(3,-1).T]
    coordinates[:count]=coordinates[:count]*unit_conversion
    f.close()
    if not nr_return == None:
        nr_return.append(nr_atoms)
        nr_return.append(len(charges)-nr_atoms)
    return coordinates,charges

#taken from http://stackoverflow.com/questions/2885190/using-pythons-configparser-to-read-a-file-without-section-name
#and modified to be less elaborate
class _SectionlessConfigParser(ConfigParser.ConfigParser):
    """
    Extends ConfigParser to allow files without sections.

    This is done by wrapping read files and prepending them with a placeholder
    section, which defaults to '__DEFAULT__'
    """

    def __init__(self, *args, **kwargs):
        ConfigParser.ConfigParser.__init__(self, *args, **kwargs)

    def readfp(self, fp, *args, **kwargs):
        import StringIO
        with open(fp) as stream:
            fakefile = StringIO.StringIO("[__DEFAULT__]\n" + stream.read())
        return ConfigParser.ConfigParser.readfp(self, fakefile, *args, **kwargs)
    
    def get_int(self,name, *args, **kwargs):
        return int(self.get("__DEFAULT__",name, *args, **kwargs))

    def get_float(self,name, *args, **kwargs):
        return float(self.get("__DEFAULT__",name, *args, **kwargs))

    def get_boolean(self,name, *args, **kwargs): 
        return self.getboolean("__DEFAULT__",name, *args, **kwargs)

    def allitems(self,name, *args, **kwargs):
        return self.items("__DEFAULT__", *args, **kwargs)

    def get_str(self,name, *args, **kwargs):
        return self.get("__DEFAULT__", name, *args, **kwargs)

def read_config_file(filename,defaults=None):
    """
    Read in a section-less config file and return an appropriate config parser object.
    Use methods "parser.get_str('name')" to get a string object corresponding to the name.
    Other functions defined: get_int, get_float, get_boolean, allitems (return all items
    as strings in a dictionary). Types appropriate to the name will be returned.

    filename: the name of the config file to read in
    """
    #info on how to do this has been taken from:
    #http://stackoverflow.com/questions/2885190/using-pythons-configparser-to-read-a-file-without-section-name
    parser = _SectionlessConfigParser(defaults=defaults)
    parser.readfp(filename)
    return parser
