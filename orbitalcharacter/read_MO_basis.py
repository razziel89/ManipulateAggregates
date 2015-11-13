"""
High-level module that allows for reading in basis set information and
transforming it in a format suitable for further processing.
"""

import copy
import math

from collection.read import read_molden

def _get_MOs_occupation(MOs,occ_func):
    MOsalpha  = [mo[3] for mo in MOs if mo[1].lower()=='alpha' and occ_func(mo[2])]
    MOsbeta   = [mo[3] for mo in MOs if mo[1].lower()=='beta'  and occ_func(mo[2])]
    OCCsalpha = [mo[2] for mo in MOs if mo[1].lower()=='alpha' and occ_func(mo[2])]
    OCCsbeta  = [mo[2] for mo in MOs if mo[1].lower()=='beta'  and occ_func(mo[2])]
    if len(MOsalpha)>0 or len(MOsbeta)>0:
        return (MOsalpha,OCCsalpha),(MOsbeta,OCCsbeta)
    else:
        raise ValueError("No molecular orbitals fulfilling the occupation function defined for either spin.")

def _gen_basis_from_GTO(Atoms,GTO,sorting='none'):
    """
    sorting: none: s,p,d,f orbitals for each atom one after the other
             spdf:  have an ss block, then an sp-block, etc.
    """
    GTO=[gto[1] for gto in GTO]
    if sorting.lower()=='spdf':
        for wanttype in ['s','p','d']:
            for (element,(x,y,z)),gto in zip(Atoms,GTO):
                #prefactor is deliberately being ignored because it is no longer
                #used by the molden file standard
                for shelltype,prefactor,nr_prim,prim in gto:
                    if shelltype == wanttype:
                        if shelltype == 's':
                            yield (x,y,z),(0,0,0),prim
                        elif shelltype == 'p':
                            yield (x,y,z),(1,0,0),prim
                            yield (x,y,z),(0,1,0),prim
                            yield (x,y,z),(0,0,1),prim
                        elif shelltype == 'd':
                            yield (x,y,z),(2,0,0),prim
                            yield (x,y,z),(0,2,0),prim
                            yield (x,y,z),(0,0,2),prim
                            yield (x,y,z),(1,1,0),prim
                            yield (x,y,z),(1,0,1),prim
                            yield (x,y,z),(0,1,1),prim
                        elif shelltype == 'f':
                            yield (x,y,z),(3,0,0),prim
                            yield (x,y,z),(0,3,0),prim
                            yield (x,y,z),(0,0,3),prim
                            yield (x,y,z),(2,1,0),prim
                            yield (x,y,z),(2,0,1),prim
                            yield (x,y,z),(1,2,0),prim
                            yield (x,y,z),(0,2,1),prim
                            yield (x,y,z),(1,0,2),prim
                            yield (x,y,z),(0,1,2),prim
                            yield (x,y,z),(1,1,1),prim
                        else:
                            raise ValueError("Only s,p,d and f-type shells are supported.")
    elif sorting.lower()=='none':
        for (element,(x,y,z)),gto in zip(Atoms,GTO):
            #prefactor is deliberately being ignored because it is no longer
            #used by the molden file standard
            for shelltype,prefactor,nr_prim,prim in gto:
                if shelltype == 's':
                    yield (x,y,z),(0,0,0),prim
                elif shelltype == 'p':
                    yield (x,y,z),(1,0,0),prim
                    yield (x,y,z),(0,1,0),prim
                    yield (x,y,z),(0,0,1),prim
                elif shelltype == 'd':
                    yield (x,y,z),(2,0,0),prim
                    yield (x,y,z),(0,2,0),prim
                    yield (x,y,z),(0,0,2),prim
                    yield (x,y,z),(1,1,0),prim
                    yield (x,y,z),(1,0,1),prim
                    yield (x,y,z),(0,1,1),prim
                elif shelltype == 'f':
                    yield (x,y,z),(3,0,0),prim
                    yield (x,y,z),(0,3,0),prim
                    yield (x,y,z),(0,0,3),prim
                    yield (x,y,z),(1,2,0),prim
                    yield (x,y,z),(2,1,0),prim
                    yield (x,y,z),(2,0,1),prim
                    yield (x,y,z),(1,0,2),prim
                    yield (x,y,z),(0,1,2),prim
                    yield (x,y,z),(0,2,1),prim
                    yield (x,y,z),(1,1,1),prim
                else:
                    raise ValueError("Only s,p,d and f-type shells are supported.")
    else:
        raise ValueError("Sorting must be None or SPD")

#all variable names that contain "alpha" or "beta" are spin-polarized values
def get_MOs_and_basis(filename, occ_func=lambda o:True, filetype="molden", spins='alpha', copy_non_present=True, msave=None, sorting='none'):
    """
    Read information about molecular orbitals and the basis from a file.
    Basis functions are contracted Cartesian Gaussian functions.

    filename: str
        The name of the file from which to read in the data
    occ_func: function, returns boolean
        A function that selects orbitals based on their occupation numbers.
        Give "lambda o: True" (without the quotes) if you want every orbital
        to be returned (default).
    filetype: str
        The type of file from which to read the data.
        Possible options: molden
    spins: str
        Declare which spins to return. Possible values: alpha, beta, both
    copy_non_present: bool
        If a requested spin is not present in the file, replace it by the
        present one. That way declaring 'both' for spins will return two
        lists even if no spin 'beta' or 'alpha' is given.
    msave: dictionary
        If a dictionary is given, add data from the internal variable m.
        This makes printing out another molden file using the same basis
    sorting: none: s,p,d,f orbitals for each atom one after the other
             spdf:  have an ss block, then an sp-block, etc. in the basis
    """

    if spins == "both":
        want_alpha = True
        want_beta  = True
    elif spins == "alpha":
        want_alpha = True 
        want_beta  = False
    elif spins == "beta":
        want_alpha = False
        want_beta  = True
    else:
        raise ValueError("Unknown spin requested.")

    if filetype=="molden":
        m = read_molden(filename,positions=True,elementnames=True,GTO=True,GTO_coefficients=True,GTO_nr_primitives=False,MO=True,MO_coefficients=True)
    else:
        raise ValueError("Unsupported filetype specified.")

    (MOsalpha,OCCsalpha),(MOsbeta,OCCsbeta) = _get_MOs_occupation(m["MO"], occ_func)
    basis = list(_gen_basis_from_GTO(m["positions"],m["GTO"]))

    #one of the two spins is defintely present, otherwise _get_MOs_occupation raised an error
    if want_alpha and len(MOsalpha)==0:
        if copy_non_present:
            MOsalpha  = copy.deepcopy(MOsbeta)
            #if one spin is duplicated, the occupations have to be halved
            OCCsbeta  = [0.5*occ for occ in OCCsbeta]
            OCCsalpha = copy.deepcopy(OCCsbeta)
    if want_beta and len(MOsbeta)==0:
        if copy_non_present:
            MOsbeta   = copy.deepcopy(MOsalpha)
            #if one spin is duplicated, the occupations have to be halved
            OCCsalpha = [0.5*occ for occ in OCCsalpha]
            OCCsbeta  = copy.deepcopy(OCCsalpha)
    
    if not want_alpha:
        del MOsalpha
        MOsalpha = None
        del OCCsalpha
        OCCsalpha = None
    if not want_beta:
        del MOsbeta
        MOsbeta = None
        del OCCsbeta
        OCCsbeta = None

    if msave is not None:
        for key in m:
            msave[key] = copy.deepcopy(m[key])

    return basis,(MOsalpha,OCCsalpha),(MOsbeta,OCCsbeta)
