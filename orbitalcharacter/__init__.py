"""
Check for orbital character.
"""

import sys
import copy

import numpy as np
import math

from collection.read import read_charges_dx as rdx
from collection.write import print_molden as pm
from collection.write import print_dx_file as pdx

from orbitalcharacter.read_MO_basis import get_MOs_and_basis
from orbitalcharacter.Smatrix import Smatrix, normalize_basis, overlap_lincomb, normalize_MO

def _expansion(MOs,Smatrix,normalize=False):
    Bsize  = len(Smatrix)
    MOsize = len(MOs)
    RHS    = [    sum((
                        coeff*Skj
                   for mo in MOs for Skj,coeff in zip(Sk,mo)
                   ))
             for Sk in Smatrix
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
    result = np.linalg.solve(Smatrix,RHS)
    if normalize:
        result /= sqrt(_overlap_lincomb(Smatrix,result))
    return result

#all variable names that contain "alpha" or "beta" are spin-polarized values
if __name__ == "__main__":
    #read in the molden files and extract spin-polarized MO information from them
    #also read in basis information
    #only occupied orbitals are being read in (occupation > 0.1)
    print >>sys.stderr,"DEBUG: started"
    m={} #save all the data of the molden file to be able to re-print everything
    basis,MOs1alpha,MOs1beta  = get_MOs_and_basis(sys.argv[1],occ_func=lambda o:o>0.1,filetype="molden",spins='both',msave=m)
    basis2,MOs2alpha,MOs2beta = get_MOs_and_basis(sys.argv[2],occ_func=lambda o:o>0.1,filetype="molden",spins='both')
    print >>sys.stderr,"DEBUG: reading molden-files and basis generation done"
    #this compares all elements in basis and basis2 recursively and only returns True
    #if all elements are the same (i.e. have the same value and are of the same shape)
    if not basis == basis2:
        raise ValueError("The bases defined in the two molden files are not the same, cannot compare.")
    del basis2
    print >>sys.stderr,"DEBUG: basis comparison done"
    Smat = Smatrix(basis)
    print >>sys.stderr,"DEBUG: built S matrix"
    Smat = normalize_basis(basis,Smat)
    print >>sys.stderr,"DEBUG: renormalized basis and corrected Smatrix"
    #print [sorted([
    #            (abs(overlap_lincomb(Smat,coeffs[i],coeffs[j])),i,j)
    #      for i in xrange(len(coeffs)) for j in xrange(len(coeffs)) if i!=j
    #      ],key=lambda o:o[0])[-10:-1] for coeffs in [MOs1alpha,MOs1beta,MOs2alpha,MOs2beta]]
    #print >>sys.stderr,"DEBUG: printed overlap"
    normalize_MO(Smat,MOs1alpha)
    normalize_MO(Smat,MOs1beta)
    normalize_MO(Smat,MOs2alpha)
    normalize_MO(Smat,MOs2beta)
    print >>sys.stderr,"DEBUG: renormalized molecular orbitals"
    #print [sorted([
    #            (abs(overlap_lincomb(Smat,coeffs[i],coeffs[j])),i,j)
    #      for i in xrange(len(coeffs)) for j in xrange(len(coeffs)) if i!=j
    #      ],key=lambda o:o[0])[-10:-1] for coeffs in [MOs1alpha,MOs1beta,MOs2alpha,MOs2beta]]
    print >>sys.stderr,"DEBUG: done"

##expand all total wave functions in terms of the basis functions
#psi1alpha = _expansion(MOs1alpha,Smatrix,normalize=False)
#psi1beta  = _expansion(MOs1beta,Smatrix,normalize=False)
#psi2alpha = _expansion(MOs2alpha,Smatrix,normalize=False)
#psi2beta  = _expansion(MOs2beta,Smatrix,normalize=False)
#
#print >>sys.stderr,"DEBUG: expansion in terms of basis functions done"
#
#MO = [[1.0,'alpha',2.0,psi1alpha],[2.0,'alpha',2.0,psi1beta],[3.0,'alpha',2.0,psi2alpha],[4.0,'alpha',2.0,psi2beta]]
#pm("test.molden",positions=m["positions"],element_names=True,GTO=m["GTO"],MO=MO)
#print >>sys.stderr,"DEBUG: printing molden file done"
#
#header={}
##invert_charge_data only True for this test case
#dxfile = rdx(sys.argv[3],invert_charge_data=True,rescale_charges=False,density=True,header_dict=header)
#grid=dxfile[0]
#print >>sys.stderr,"DEBUG: reading example dx-file done"
#
#data = _prepare_grid_calculation(grid,basis)
#print >>sys.stderr,"DEBUG: prepared grid calculation"
#
#rho1alpha,rho1beta,rho2alpha,rho2beta = _density_on_grid([psi1alpha,psi1beta,psi2alpha,psi2beta],*data)
#print >>sys.stderr,"DEBUG: generated density on grid"
#
#pdx("rho1alpha.dx",header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],rho1alpha)
#pdx("rho1beta.dx",header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],rho1beta)
#pdx("rho2alpha.dx",header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],rho2alpha)
#pdx("rho2bega.dx",header["counts_xyz"],header["org_xyz"],header["delta_x"],header["delta_y"],header["delta_z"],rho2beta)
#print >>sys.stderr,"DEBUG: wrote dx-files"
