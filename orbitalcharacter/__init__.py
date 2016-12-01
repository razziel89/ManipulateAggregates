"""This submodule agregates functions that can be applied to QM orbitals.

@package ManipulateAggregates.orbitalcharacter
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
import sys
import copy

import logging
logger = logging.getLogger(__name__)
try:
    import numpy as numpy
except ImportError:
    logger.warning("Could not import numpy")

try:
    from . import density_on_grid as density_on_grid
except ImportError:
    logger.warning("Could not import .density_on_grid")
try:
    from . import read_MO_basis   as read_MO_basis
except ImportError:
    logger.warning("Could not import .read_MO_basis")
try:
    from . import Smatrix         as Smatrix
except ImportError:
    logger.warning("Could not import .Smatrix")
    
try:
    from ..collection.write import print_dx_file
except ImportError:
    logger.warning("Could not import ..collection.write.print_dx_file")
try:
    from ..collection.read import read_dx
except ImportError:
    logger.warning("Could not import ..collection.read.read_dx")
try:
    from ..collection.read import read_molden
except ImportError:
    logger.warning("Could not import ..collection.read.read_molden")
    
try:
    from ..orbitalcharacter.read_MO_basis import get_MOs_and_basis
except ImportError:
    logger.warning("Could not import ..orbitalcharacter.read_MO_basis.get_MOs_and_basis")
try:
    from ..orbitalcharacter.Smatrix import Smatrix, normalize_basis, overlap_lincomb, normalize_MOs
except ImportError:
    logger.warning("Could not import Smatrix, normalize_basis, overlap_lincomb or normalize_MOs from ..orbitalcharacter.Smatrix")
    
try:
    from FireDeamon import ElectronDensityPy, InitializeGridCalculationOrbitalsPy
except ImportError:
    logger.warning("Could not import ElectronDensityPy or InitializeGridCalculationOrbitalsPy from FireDeamon")
try:
    from FireDeamon import InitializeGridCalculationOrbitalsPy, ElectrostaticPotentialOrbitalsPy, ElectrostaticPotentialPy
except ImportError:
    logger.warning("Could not import InitializeGridCalculationOrbitalsPy, ElectrostaticPotentialOrbitalsPy or ElectrostaticPotentialPy from FireDeamon")

## all orbitals with a higher occupation than this are considered to be occupied
MINOCC = 0.0000001
## convertion factor from Bohr (atomic units) to Angstroms
BOHRTOANG = 0.529177249

def expand_total_wavefunction(MOs,Smat,normalize=False):
    """Compute the total wavefunction, i.e., the sum over all molecular orbitals.

    Args:
        MOs: (list of lists of floats) coefficients describing the molecular orbitals
        Smat: (square matrix of floats) matrix describing the overlap
            between the shels of the basis used to obtain the molecular orbitals

    Kwargs:
        normalize: (bool) whether or not to make it so that the overlap of the
            returned total wavefunction with itself shall be 1

    Returns:
        a list of floats describing the total wavefunction in the given basis
    """
    Bsize  = len(Smat)
    MOsize = len(MOs)
    RHS    = [    sum((
                        coeff*Skj
                   for mo in MOs for Skj,coeff in zip(Sk,mo)
                   ))
             for Sk in Smat
             ]
    #I want to have: |P> = SUM (d_k*|phi_k>) for i in inteval [0,N]
    #with: P: total wavefunction, d_k: linear combination coefficient
    #      S_ik: element of the overlap matrix i.e. <phi_i|phi_k> (symmetric)
    #      N: nr. basis functions
    #      |phi_k>: k-th basis function
    #Hence, <phi_i|P> = SUM (d_k <phi_i|phi_k>) for k in interval [0,N]
    #                 = SUM (d_k S_ik) for k in interval [0,N]
    #With v = [<phi_1|P>,<phi_2|P>,...,<phi_N|P>] and d = [d_1,d_2,...,d_N] (Python lists as vectors)
    #     follows v = S dot d (dot: matrix product)
    #     and hence the equation S*d=v has to be solved for v
    #solves S*d=RHS where d is the vector containing the coefficients of interest
    result = numpy.linalg.solve(Smat,RHS)
    if normalize:
        result /= sqrt(overlap_lincomb(Smat,result))
    return result

def read_molden_orbitals_corrected(filename):
    """Read in a molden-file and apply corrections.

    The applied corrections make sure that each shell in the given basis is
    normalized and that all molecular orbitals are normalized.

    Args:
        filename: (string) the name of the molden-file ro be read. Not a path.

    Returns:
        basis,Smat,(MOsalpha,MOsbeta),(OCCsalpha,OCCsbeta). The value for basis
        (a list of [A,L,Prim]) is what is described
        in ManipulateAggregates.orbitalcharacter.density_on_grid.basis. Smat is
        described in ManipulateAggregates.orbitalcharacter.expand_total_wavefunction.
        MOsalpha and MOsbeta are lists of floats of the molecular orbital
        coefficients for alpha and beta spins, respectively. OCCsalpha and
        OCCsbeta are lists of floats of the occupations of the molecular
        orbitals for alpha and beta spins, respectively.
        
    """
    print >>sys.stderr,"DEBUG: reading in basis from molden file and applying corrections for limited precision read"
    #read in the molden file and extract spin-polarized MO information from it
    #also read in basis information
    occ_func=lambda o: o>MINOCC
    basis,(allMOsalpha,allOCCsalpha),(allMOsbeta,allOCCsbeta),(IdxHOMOalpha,IdxHOMObeta)  = get_MOs_and_basis(
            filename,filetype="molden",spins='both',
            alpha_high_energy=True,occ_func=occ_func)
    #determine occupied orbitals
    MOsalpha  = allMOsalpha [:IdxHOMOalpha+1]
    MOsbeta   = allMOsbeta  [:IdxHOMObeta +1]
    OCCsalpha = allOCCsalpha[:IdxHOMOalpha+1]
    OCCsbeta  = allOCCsbeta [:IdxHOMObeta +1]
    Smat = Smatrix(basis)
    Smat = normalize_basis(basis,Smat) #after this, basis will be normalized
    copy_beta = (MOsalpha==MOsbeta)
    normalize_MOs(Smat,MOsalpha,occupations=OCCsalpha)
    if copy_beta:
        MOsbeta = copy.deepcopy(MOsalpha)
    else:
        normalize_MOs(Smat,MOsbeta,occupations=OCCsbeta)
    return basis,Smat,(MOsalpha,MOsbeta),(OCCsalpha,OCCsbeta)

def edensity(filename,header=None,dir="",progress=False,save_mos=(0,0),points=80,cutoff=7.0,gzipped=True,special_orbitals=True,outfile="rho.dx"):
    """Compute the total (and optionally partial) electron density on a regular grid.

    If @a header is None, a suitable header will be autogenerated.

    Args:
        filename: (string) the name of the input molden-file

    Kwargs:
        header: (dictionary) has to have the keys "counts_xyz", "org_xyz",
            "delta_x", "delta_y", "delta_z" all present. There meanings are what is
            mentioned in ManipulateAggregates.collection.write.print_dx_file. All
            values have to be in Angstroms.
        dir: (string) directory name to be prefixed to all output files.
        progress: (bool) whether or not progress shall be reported during the computation
        save_mos: (tuple of 2 int) a range of which orbitals shall be printed
            separately to files.  Counting starts at 1. If 0 is provided, it
            means "last occupied orbital" (only if the other value is not 0).
        points: (int) how many points shall be used in each Cartesian direction
            if the grid has to be automatically generated.
        cutoff: (float, in atomic units!!!) if a gridpoint and the center of a
            basis function are farther apart from each other than this value,
            the density will not be evaluated.
        gzipped: (bool) whether or not the dx-files shall be written in gzip
            compresed format.
        special_orbitals: (bool) whether or not the density of some "special"
            orbitals shall always be output. Special orbitals are the frontier
            orbitals, i.e., HOMO and LUMO.
        outfile: (string) to which file the total density shall be written. Not
            a path.
    """
    if len(dir)>0:
        if not dir.endswith("/"):
            dir+="/"
    #read in the molden file and extract spin-polarized MO information from it
    #also read in basis information
    occ_func=lambda o: o>MINOCC
    print >>sys.stderr,"DEBUG: started comptation of electronic density"
    basis,(allMOsalpha,allOCCsalpha),(allMOsbeta,allOCCsbeta),(IdxHOMOalpha,IdxHOMObeta)  = get_MOs_and_basis(
            filename,filetype="molden",spins='both',
            alpha_high_energy=True,occ_func=occ_func)
    #determine occupied orbitals
    MOsalpha  = allMOsalpha [:IdxHOMOalpha+1]
    MOsbeta   = allMOsbeta  [:IdxHOMObeta +1]
    OCCsalpha = allOCCsalpha[:IdxHOMOalpha+1]
    OCCsbeta  = allOCCsbeta [:IdxHOMObeta +1]
    #determine LUMOs (index starts at 0, so the length is always maxindex+1)
    LUMOalpha = allMOsalpha[IdxHOMOalpha+1]
    LUMObeta  = allMOsalpha[IdxHOMOalpha+1]
    print >>sys.stderr,"DEBUG: reading molden-files and basis generation done"
    if header is None:
        #read_molden returns coordinates in Bohr but I want Angstroms at this position
        coordinates = numpy.array([c[1] for c in read_molden(filename,positions=True,GTO=False,MO=False)["positions"]])*BOHRTOANG
        min_corner = numpy.amin(coordinates,axis=0)-10.0
        max_corner = numpy.amax(coordinates,axis=0)+10.0
        counts_xyz = numpy.array([points,points,points])
        org_xyz   = min_corner
        #grid creation copied from energyscan.scan but slightly altered
        space = [numpy.linspace(s,e,num=c,dtype=float)
                    for s,e,c
                    in zip(min_corner,max_corner,counts_xyz)
               ]
        #just take the difference between the first elements in every direction to get the stepsize
        delta_x = numpy.array([space[0][1] - space[0][0], 0.0, 0.0])
        delta_y = numpy.array([0.0, space[1][1] - space[1][0], 0.0])
        delta_z = numpy.array([0.0, 0.0, space[2][1] - space[2][0]])
        a1,a2,a3  = numpy.array(numpy.meshgrid(*space,indexing="ij"))
        a1.shape  = (-1,1)
        a2.shape  = (-1,1)
        a3.shape  = (-1,1)
        grid      = numpy.concatenate((a1,a2,a3),axis=1)
        print >>sys.stderr,"DEBUG: autogenerated header and grid for dx-files"
    else:
        counts_xyz = numpy.array(header["counts_xyz"])
        org_xyz    = numpy.array(header["org_xyz"])
        delta_x    = numpy.array(header["delta_x"])
        delta_y    = numpy.array(header["delta_y"])
        delta_z    = numpy.array(header["delta_z"])
        print >>sys.stderr,"DEBUG: done reading data from header"
        #grid creation copied from energyscan.py but slightly altered
        space = [numpy.linspace(s,e,num=c,dtype=float)
                    for s,e,c
                    in zip(org_xyz,org_xyz+counts_xyz[0]*delta_x+counts_xyz[1]*delta_y+counts_xyz[2]*delta_z,counts_xyz)
               ]
        #just take the difference between the first elements in every direction to get the stepsize
        a1,a2,a3  = numpy.array(numpy.meshgrid(*space,indexing="ij"))
        a1.shape  = (-1,1)
        a2.shape  = (-1,1)
        a3.shape  = (-1,1)
        grid      = numpy.concatenate((a1,a2,a3),axis=1)
        print >>sys.stderr,"DEBUG: generated grid for dx-files from header"
    Smat = Smatrix(basis)
    print >>sys.stderr,"DEBUG: built S matrix"
    Smat = normalize_basis(basis,Smat) #after this, basis will be normalized
    print >>sys.stderr,"DEBUG: renormalized basis and corrected Smatrix"
    normalize_MOs(Smat,MOsalpha,occupations=OCCsalpha)
    normalize_MOs(Smat,MOsbeta,occupations=OCCsbeta)
    print >>sys.stderr,"DEBUG: renormalized molecular orbitals"
    #expand all total wave functions in terms of the basis functions
    psialpha = expand_total_wavefunction(MOsalpha,Smat)
    psibeta  = expand_total_wavefunction(MOsbeta,Smat)
    print >>sys.stderr,"DEBUG: expansion in terms of basis functions done"
    nr_electrons_alpha = int(round(overlap_lincomb(Smat,psialpha)))
    nr_electrons_beta  = int(round(overlap_lincomb(Smat,psibeta)))
    nr_electrons       = nr_electrons_alpha + nr_electrons_beta
    #the grid has been read in in Angstroms so it has to be converted to bohrs
    #data = prepare_grid_calculation(grid,basis,scale=0.52918) 
    data = InitializeGridCalculationOrbitalsPy(grid,basis,scale=BOHRTOANG) 
    print >>sys.stderr,"DEBUG: prepared grid calculation"
    async = progress
    save_all_mos = not(save_mos == (0,0))
    if save_all_mos:
        if save_mos[0] == 0:
            save_mos = (max((IdxHOMOalpha,IdxHOMObeta))+1,save_mos[1])
        if save_mos[1] == 0:
            save_mos = (save_mos[0],max((IdxHOMOalpha,IdxHOMObeta))+1)
        if save_mos[0] > save_mos[1]:
            save_mos = (save_mos[1],save_mos[0])
        tot_dens = numpy.zeros((len(grid),),dtype=float)
        if MOsalpha == MOsbeta:
            mocount=1
            for mo,occ in zip(MOsalpha,OCCsalpha):
                tempdens = numpy.array(ElectronDensityPy([mo],data,occupations=[occ],cutoff=cutoff,prog_report=async))
                if mocount >= save_mos[0] and mocount <= save_mos[1]:
                    print_dx_file(dir+"MO"+str(mocount)+"alpha.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,tempdens,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
                    print_dx_file(dir+"MO"+str(mocount)+"beta.dx", counts_xyz,org_xyz,delta_x,delta_y,delta_z,tempdens,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
                tot_dens += 2*tempdens
                mocount += 1
        else:
            mocount=1
            for mo,occ in zip(MOsalpha,OCCsalpha):
                tempdens = numpy.array(ElectronDensityPy([mo],data,occupations=[occ],cutoff=cutoff,prog_report=async))
                if mocount >= save_mos[0] and mocount <= save_mos[1]:
                    print_dx_file(dir+"MO"+str(mocount)+"alpha.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,tempdens,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
                tot_dens += tempdens
                mocount += 1
            mocount=1
            for mo,occ in zip(MOsbeta,OCCsbeta):
                tempdens = numpy.array(ElectronDensityPy([mo],data,occupations=[occ],cutoff=cutoff,prog_report=async))
                if mocount >= save_mos[0] and mocount <= save_mos[1]:
                    print_dx_file(dir+"MO"+str(mocount)+"beta.dx", counts_xyz,org_xyz,delta_x,delta_y,delta_z,tempdens,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
                tot_dens += tempdens
                mocount += 1
    else:
        if MOsalpha == MOsbeta:
            tot_dens = 2.0*numpy.array(ElectronDensityPy(MOsalpha,data,occupations=OCCsalpha,cutoff=cutoff,prog_report=async))
        else:
            tot_dens = numpy.array(ElectronDensityPy(MOsalpha+MOsbeta,data,occupations=OCCsalpha+OCCsbeta,cutoff=cutoff,prog_report=async))
    print >>sys.stderr,"DEBUG: generated total density on grid"
    print_dx_file(dir+outfile,counts_xyz,org_xyz,delta_x,delta_y,delta_z,tot_dens,comment="Nr. Electrons: %d"%(nr_electrons),gzipped=gzipped)
    print >>sys.stderr,"DEBUG: wrote dx-file for total density"
    if special_orbitals:
        dens_homo_alpha = numpy.array(ElectronDensityPy([MOsalpha[-1]],data,occupations=[OCCsalpha[-1]],cutoff=cutoff,prog_report=async))
        dens_homo_beta  = numpy.array(ElectronDensityPy([MOsbeta[-1]],data,occupations=[OCCsbeta[-1]],cutoff=cutoff,prog_report=async))
        print >>sys.stderr,"DEBUG: computed HOMO densities (both spins)"
        print_dx_file(dir+"HOMO1.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,dens_homo_alpha,comment="Nr. Electrons: %d"%(OCCsalpha[-1]),gzipped=gzipped)
        print_dx_file(dir+"HOMO2.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,dens_homo_beta,comment="Nr. Electrons: %d"%(OCCsbeta[-1]),gzipped=gzipped)
        print >>sys.stderr,"DEBUG: wrote dx-files for HOMO densities (both spins)"
        dens_lumo_alpha = numpy.array(ElectronDensityPy([LUMOalpha],data,occupations=[1.0],cutoff=cutoff,prog_report=async))
        dens_lumo_beta  = numpy.array(ElectronDensityPy([LUMOalpha],data,occupations=[1.0],cutoff=cutoff,prog_report=async))
        print >>sys.stderr,"DEBUG: computed LUMO densities (both spins)"
        #occupation is not defined for the LUMO so normalization does not work (normalize to 1)
        print_dx_file(dir+"LUMO1.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,dens_lumo_alpha,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
        print_dx_file(dir+"LUMO2.dx",counts_xyz,org_xyz,delta_x,delta_y,delta_z,dens_lumo_beta,comment="Nr. Electrons: %d"%(1),gzipped=gzipped)
        print >>sys.stderr,"DEBUG: wrote dx-files for LUMO densities (both spins)"
    print >>sys.stderr,"DEBUG: done computation of electronci density"

def _correlation(array1,array2):
    """Return the correlation between two NumPy arrays."""
    temparray1 = array1-numpy.mean(array1)
    temparray2 = array2-numpy.mean(array2)
    temparray1 = temparray1/numpy.sqrt(numpy.dot(temparray1,temparray1))
    temparray2 = temparray2/numpy.sqrt(numpy.dot(temparray2,temparray2))
    return numpy.dot(temparray1,temparray2)

def electrostatic_potential(filename,header=None,dir="",progress=False,points=80,
        ext_grid=None,at_coordinates=None,charges=None,outfile="potential.dx"):
    """Compute the total electrostatic potential on a regular grid.

    If @a header is None, a suitable header will be autogenerated. The
    computation only makes sense if also the charges of the molecule's atoms
    are considered. These are given via @a at_coordinates and @a charges.

    Args:
        filename: (string) the name of the input molden-file. Only the [GTO]
            and [MO] sections are required.

    Kwargs:
        header: (dictionary) has to have the keys "counts_xyz", "org_xyz",
            "delta_x", "delta_y", "delta_z" all present. There meanings are what is
            mentioned in ManipulateAggregates.collection.write.print_dx_file. All
            values have to be in Angstroms.
        dir: (string) directory name to be prefixed to all output files.
        progress: (bool) whether or not progress shall be reported during the computation
        points: (int) how many points shall be used in each Cartesian direction
            if the grid has to be automatically generated.
        ext_grid: (list of lists of 3 floats, length divisible by 3) the
            externally given grid on which the electrostatic potential shall
            be computed.
        at_coordinates: (list of lists of 3 floats) each sublist contains the
            Cartesian coordinates of the atoms.
        charges: (list of floats) charges associated with the atoms
        outfile: (string) to which file the total potential shall be written.
            Not a path.
    """
    if ext_grid is not None:
        raise Exception("EXTERNAL GRID NOT YET IMPLEMENTED")
    if len(dir)>0:
        if not dir.endswith("/"):
            dir+="/"
    #read in the molden file and extract spin-polarized MO information from it
    #also read in basis information
    MINOCC = 0.1 #all orbitals with a higher occupation than this are considered to be occupied
    occ_func=lambda o: o>MINOCC
    print >>sys.stderr,"DEBUG: started computation of electrostatic potential"
    basis,(allMOsalpha,allOCCsalpha),(allMOsbeta,allOCCsbeta),(IdxHOMOalpha,IdxHOMObeta)  = get_MOs_and_basis(
            filename,filetype="molden",spins='both',
            alpha_high_energy=True,occ_func=occ_func)
    #determine occupied orbitals
    MOsalpha  = allMOsalpha [:IdxHOMOalpha+1]
    MOsbeta   = allMOsbeta  [:IdxHOMObeta +1]
    OCCsalpha = allOCCsalpha[:IdxHOMOalpha+1]
    OCCsbeta  = allOCCsbeta [:IdxHOMObeta +1]
    #determine LUMOs (index starts at 0, so the length is always maxindex+1)
    LUMOalpha = allMOsalpha[IdxHOMOalpha+1]
    LUMObeta  = allMOsalpha[IdxHOMOalpha+1]
    print >>sys.stderr,"DEBUG: reading molden-files and basis generation done"
    if header is None:
        #read_molden returns coordinates in Bohr but I want Angstroms at this position
        coordinates = numpy.array([c[1] for c in read_molden(filename,positions=True,GTO=False,MO=False)["positions"]])*BOHRTOANG
        min_corner = (numpy.amin(coordinates,axis=0)/BOHRTOANG-10.0)*BOHRTOANG
        max_corner = (numpy.amax(coordinates,axis=0)/BOHRTOANG+10.0)*BOHRTOANG
        #min_corner = numpy.array([0.5,1,0.5],dtype=float)
        #max_corner = numpy.array([0.5,4,0.5],dtype=float)
        #counts_xyz = numpy.array([1,4,1])
        #min_corner = numpy.amin(coordinates,axis=0)-10.0
        #max_corner = numpy.amax(coordinates,axis=0)+10.0
        counts_xyz = numpy.array([points,points,points])
        org_xyz    = min_corner
        #grid creation copied from energyscan.scan but slightly altered
        space = [numpy.linspace(s,e,num=c,dtype=float)
                    for s,e,c
                    in zip(min_corner,max_corner,counts_xyz)
               ]
        #just take the difference between the first elements in every direction to get the stepsize
        delta_x = numpy.array([space[0][1] - space[0][0], 0.0, 0.0])
        delta_y = numpy.array([0.0, space[1][1] - space[1][0], 0.0])
        delta_z = numpy.array([0.0, 0.0, space[2][1] - space[2][0]])
        #delta_x = numpy.array([0.0, 0.0, 0.0])
        #delta_y = numpy.array([0.0, space[1][1] - space[1][0], 0.0])
        #delta_z = numpy.array([0.0, 0.0, 0.0])
        a1,a2,a3  = numpy.array(numpy.meshgrid(*space,indexing="ij"))
        a1.shape  = (-1,1)
        a2.shape  = (-1,1)
        a3.shape  = (-1,1)
        grid      = numpy.concatenate((a1,a2,a3),axis=1)
        #print grid
        print >>sys.stderr,"DEBUG: autogenerated header and grid for dx-files"
    else:
        counts_xyz = numpy.array(header["counts_xyz"])
        org_xyz    = numpy.array(header["org_xyz"])
        delta_x    = numpy.array(header["delta_x"])
        delta_y    = numpy.array(header["delta_y"])
        delta_z    = numpy.array(header["delta_z"])
        print >>sys.stderr,"DEBUG: done reading data from header"
        #grid creation copied from energyscan.py but slightly altered
        space = [numpy.linspace(s,e,num=c,dtype=float)
                    for s,e,c
                    in zip(org_xyz,org_xyz+counts_xyz[0]*delta_x+counts_xyz[1]*delta_y+counts_xyz[2]*delta_z,counts_xyz)
               ]
        #just take the difference between the first elements in every direction to get the stepsize
        a1,a2,a3  = numpy.array(numpy.meshgrid(*space,indexing="ij"))
        a1.shape  = (-1,1)
        a2.shape  = (-1,1)
        a3.shape  = (-1,1)
        grid      = numpy.concatenate((a1,a2,a3),axis=1)
        print >>sys.stderr,"DEBUG: generated grid for dx-files from header"
    Smat = Smatrix(basis)
    print >>sys.stderr,"DEBUG: built S matrix"
    Smat = normalize_basis(basis,Smat) #after this, basis will be normalized
    print >>sys.stderr,"DEBUG: renormalized basis and corrected Smatrix"
    normalize_MOs(Smat,MOsalpha,occupations=OCCsalpha)
    normalize_MOs(Smat,MOsbeta,occupations=OCCsbeta)
    print >>sys.stderr,"DEBUG: renormalized molecular orbitals"
    #expand all total wave functions in terms of the basis functions
    psialpha = expand_total_wavefunction(MOsalpha,Smat)
    psibeta  = expand_total_wavefunction(MOsbeta,Smat)
    print >>sys.stderr,"DEBUG: expansion in terms of basis functions done"
    nr_electrons_alpha = int(round(overlap_lincomb(Smat,psialpha)))
    nr_electrons_beta  = int(round(overlap_lincomb(Smat,psibeta)))
    nr_electrons       = nr_electrons_alpha + nr_electrons_beta
    print >>sys.stderr,"DEBUG: computed number of electrons"
    #the grid has been read in in Angstroms so it has to be converted to bohrs
    data = InitializeGridCalculationOrbitalsPy(grid,basis,scale=BOHRTOANG)
    print >>sys.stderr,"DEBUG: prepared grid calculation"
    if MOsalpha == MOsbeta:
        potential = -numpy.array(ElectrostaticPotentialOrbitalsPy(MOsalpha,Smat,[2*o for o in OCCsalpha],data,prog_report=progress))
    else:
        potential = -numpy.array(ElectrostaticPotentialOrbitalsPy(MOsalpha+MOsbeta,Smat,OCCsalpha+OCCsbeta,data,prog_report=progress))
    if at_coordinates is not None and charges is not None:
        #also consider core charges
        pospotential = numpy.array(ElectrostaticPotentialPy(grid/BOHRTOANG, charges, [[xyz/BOHRTOANG for xyz in a] for a in at_coordinates]))
        #print_dx_file(dir+"neg_"+outfile,counts_xyz,org_xyz,delta_x,delta_y,delta_z,potential,gzipped=False)
        #print_dx_file(dir+"pos_"+outfile,counts_xyz,org_xyz,delta_x,delta_y,delta_z,pospotential,gzipped=False)
        #print potential
        potential += pospotential
    print >>sys.stderr,"DEBUG: generated potential on grid"
    print_dx_file(dir+outfile,counts_xyz,org_xyz,delta_x,delta_y,delta_z,potential,gzipped=False)
    print >>sys.stderr,"DEBUG: wrote dx-file for potential"
    print >>sys.stderr,"DEBUG: done computation of electrostatic potential"

def _similarity(diffdens,MOdens,type=0,name=False):
    """Compute the similarity between a difference density
    and the density of a molecular orbital.
    """
    if type == 0:
        #This should be as close to 1 as possible
        #but I guess the other two are a better
        #meassure
        result = _correlation(diffdens,MOdens)
        rname="Correlation"
    elif type == 1:
        #This should be as close to 1 as possible
        #which means that a lot of density is being 
        #taken from the MO
        temparray1 = numpy.copy(diffdens)
        temparray2 = MOdens
        temparray1[temparray1<0.0] = 0.0
        result = _correlation(temparray1,temparray2)
        rname="Correlation Positive"
    elif type == 2:
        #This should be as close to 0 as possible
        #which means that no density is being transferred
        #into the MO
        temparray1 = numpy.copy(diffdens)
        temparray2 = MOdens
        temparray1[temparray1>0.0] = 0.0
        temparray1 *= -1
        result = _correlation(temparray1,temparray2)
        rname="Correlation Negative"
    else:
        raise ValueError("Wrong type of similarity measure given.")
    if name:
        return result,rname
    else:
        return result

def density_overlap(density_1,density_2):
    """Just compute and print the correlation between two densities given by their dx-files.

    This obviously only makes sense of the two dx-files define the exact same
    grids. However, only agreement between the number of points is checked.

    Bug:
        If both dx-files define unequal grids that have the same number of
        points, the correlation is still computed but the results are not the
        actual correlations.

    Raises:
        ValueError.

    Args:
        density_1: (string) name of the dx-file that contains the first density
        density_2: (string) name of the dx-file that contains the second density
    """
    print >>sys.stderr,"DEBUG: started computation of density correlation"
    #read in all dx files
    data1  = numpy.array(read_dx(density_1,density=True,silent=True,grid=False,gzipped=True)["data"])
    data2  = numpy.array(read_dx(density_2,density=True,silent=True,grid=False,gzipped=True)["data"])
    print >>sys.stderr,"DEBUG: reading dx-files done"
    if data1.shape!=data2.shape:
        raise ValueError("Both dx files contain grids with a different number of points.")
    corrvalue = _similarity(data1,data2,0,name=False)
    print >>sys.stderr,"DEBUG: computed overlap between densities"
    print "Overlap: %8.4e"%(corrvalue)
    print >>sys.stderr,"DEBUG: done computation of density correlation"

def difference_density(density_1,density_2,dxdiffdens,compress=False,factor=1.0):
    """Just compute the difference between two densities given by their dx-files.

    Args:
        density_1: (string) name of the dx-file that contains the first density
        density_2: (string) name of the dx-file that contains the second density
        dxdiffdens: (string) name of the dx-file that shall contain the difference density.

    Kwargs:
        compress: (bool) whether or not the difference density shall be written
            in gzipped format or not.
        factor: (float) this factor is multiplied with the second density
            before computing the difference
    """
    print >>sys.stderr,"DEBUG: started computation of difference density"
    apply_func_density(density_1,density_2,dxdiffdens,compress=compress,func=lambda d1,d2:d1-factor*d2,verbose=True)
    print >>sys.stderr,"DEBUG: done computation of difference density"

global DEFAULT_FUNC
## default function for ManipulateAggregates.orbitalcharacter.apply_func_density
DEFAULT_FUNC = lambda d1,d2:d1
def apply_func_density(density_1,density_2,outdens,compress=False,func=DEFAULT_FUNC,verbose=False):
    """Apply a given function to two densities and write result to file.

    If @a density_2 is None, apply the given function only to the first
    density. If no function is provided, only the first density is written out.

    Bug:
        If both dx-files define unequal grids that have the same number of
        points, the correlation is still computed but the results are not the
        actual correlations.

    Raises:
        ValueError.

    Args:
        density_1: (string) name of the dx-file that contains the first density
        density_2: (string) name of the dx-file that contains the second density
        outdens: (string) name of the dx-file that shall contain the output density.

    Kwargs:
        compress: (bool) whether or not the difference density shall be written
            in gzipped format or not.
        func: (function of 2 variables applicable to numpy arrays) How to obtain
            the new density when given the old ones.
        verbose: (bool) give progress updates
    """
    header = {}
    #read in all dx files
    data1  = numpy.array(read_dx(density_1,density=True,silent=True,grid=False,header_dict=header,gzipped=True)["data"])
    if density_2 is not None:
        data2  = numpy.array(read_dx(density_2,density=True,silent=True,grid=False,                   gzipped=True)["data"])
    if verbose:
        print >>sys.stderr,"DEBUG: reading dx-files done"
    if density_2 is not None:
        if data1.shape!=data2.shape:
            raise ValueError("Both dx files contain grids with a different number of points.")
        newdens = func(data1,data2)
    else:
        newdens = func(data1)
    if verbose:
        print >>sys.stderr,"DEBUG: computed new density, sum: %8.4f, sum over abs: %8.4f"%(numpy.sum(diffdens),numpy.sum(numpy.fabs(diffdens)))
    print_dx_file(outdens,header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],newdens,gzipped=compress)
    if verbose:
        print >>sys.stderr,"DEBUG: wrote new density"

def _postprocess_multiple(total_1,total_2,MOalpha_1,MObeta_1,MOalpha_2,MObeta_2,dir="",type="kation"):
    """
    After creating all dx-files for each sub-calculation, use this function on
    all important dx-files to aggregate the data.

    total_1, total_2: str
        Names of the dx-files that contain the total density.
    MOalpha_1, MObeta_1: str
        Names of the dx-files that contain the HOMO density of
        alpha and beta spin, first molecule.
    MOalpha_2, MObeta_2: str
        Names of the dx-files that contain the HOMO density of
        alpha and beta spin, second molecule.
    dir: str
        Directory name to be prefixed to all output files.
    type: str
        Accepted strings are 'kation' and 'anion' depending on
        whether the neutral molecule shall be compared to the
        kation or anion.
    """
    print >>sys.stderr,"DEBUG: started postprocessing"
    if len(dir)>0:
        if not dir.endswith("/"):
            dir+="/"
    #read in all dx files
    #they have been normalized to the number of electrons
    header={}
    data1  = numpy.array(read_dx(total_1,density=True,silent=True,grid=False,header_dict=header,gzipped=True)["data"])
    data2  = numpy.array(read_dx(total_2,density=True,silent=True,grid=False                   ,gzipped=True)["data"])
    print >>sys.stderr,"DEBUG: reading dx-files done"
    nr_electrons_1 = int(round(numpy.sum(data1)))
    nr_electrons_2 = int(round(numpy.sum(data2)))
    if type=='kation':
        kation=True
        diff_to_neut=+1
        prefix="kat"
    elif type=='anion':
        kation=False
        diff_to_neut=-1
        prefix="an"
    else:
        raise ValueError("Wrong type of molecule comparison given.")
    if nr_electrons_1 == nr_electrons_2+diff_to_neut:
        neut_total  = data1
        ion_total   = data2
        if type=='kation':
            MOalpha = numpy.array(read_dx(MOalpha_1,density=True,silent=True,grid=False,gzipped=True)["data"])
            MObeta  = numpy.array(read_dx(MObeta_1, density=True,silent=True,grid=False,gzipped=True)["data"])
        else:
            MOalpha = numpy.array(read_dx(MOalpha_2,density=True,silent=True,grid=False,gzipped=True)["data"])
            MObeta  = numpy.array(read_dx(MObeta_2, density=True,silent=True,grid=False,gzipped=True)["data"])
        nr_electrons_neut = nr_electrons_1
    elif nr_electrons_1 == nr_electrons_2-diff_to_neut:
        neut_total  = data2
        ion_total   = data1
        if type=='kation':
            MOalpha = numpy.array(read_dx(MOalpha_2,density=True,silent=True,grid=False,gzipped=True)["data"])
            MObeta  = numpy.array(read_dx(MObeta_2, density=True,silent=True,grid=False,gzipped=True)["data"])
        else:
            MOalpha = numpy.array(read_dx(MOalpha_1,density=True,silent=True,grid=False,gzipped=True)["data"])
            MObeta  = numpy.array(read_dx(MObeta_1, density=True,silent=True,grid=False,gzipped=True)["data"])
        nr_electrons_neut = nr_electrons_2
    else:
        raise ValueError("Both dx files contain data about molecules that do not differ in exactly one electron.")
    print >>sys.stderr,"DEBUG: determined %sion and neutral molecule"%(prefix)
    if neut_total.shape!=ion_total.shape:
        raise ValueError("Both dx files contain grids with a different number of points.")
    diffdens = (neut_total - ion_total)*diff_to_neut
    print >>sys.stderr,"DEBUG: computed difference density, sum: %8.4f, sum over abs: %8.4f"%(numpy.sum(diffdens),numpy.sum(numpy.fabs(diffdens)))
    print_dx_file(dir+"diff_%sion.dx"%(prefix),header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],diffdens,gzipped=True)
    print >>sys.stderr,"DEBUG: wrote dx-file for difference density"
    otypes = 3
    overlap_alpha = tuple(_similarity(diffdens,MOalpha,t,name=True) for t in xrange(otypes))
    overlap_beta  = tuple(_similarity(diffdens,MObeta,t,name=True) for t in xrange(otypes))
    print >>sys.stderr,"DEBUG: computed overlap between difference density and HOMO densities (both spins)"
    for (a,na),(b,nb) in zip(overlap_alpha,overlap_beta):
        print "Type %15s: Overlap alpha/beta: %8.4e /%8.4e"%(na,a,b)
    print >>sys.stderr,"DEBUG: done postprocessing"
