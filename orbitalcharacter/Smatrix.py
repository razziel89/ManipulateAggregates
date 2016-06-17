"""
Taken from http://www.mathematica-journal.com/2012/02/evaluation-of-gaussian-molecular-integrals/
and altered a bit before transforming to Python code. All thanks goes to the authors.
"""

from collections import Sequence
from itertools import chain, count
import math

_exp   = math.exp
_sqrt  = math.sqrt

try:
    from FireDeamon import normalization_coefficient as NormCoeffPy
    from FireDeamon import Sxyz as SxyzPy
    useFD = True
except ImportError:
    useFD = False
    from scipy.special import binom as binomial

def _depth(seq):
    """
    Calculate the depth of a list. Taken from
    http://stackoverflow.com/questions/6039103/counting-deepness-or-the-deepest-level-a-nested-list-goes-to
    """
    seq = iter(seq)
    try:
        for level in count():
            seq = chain([next(seq)], seq)
            seq = chain.from_iterable(s for s in seq if isinstance(s, Sequence))
    except StopIteration:
        return level

def dfac(i):
    """
    Calculate the double factorial (not found in module "math")

    i: int
        Integer number of which to calculate the double factorial
    """
    return reduce(int.__mul__,xrange(i,0,-2),1)

def normalization_coeff(alpha,l,m,n):
    """
    Calculate the normalization coefficient of a 3D Cartesian Gaussian function
    times pi^0.75.

    alpha: float
        Exponential factor of the Gaussian
    l,m,n: each one int
        The list [l,m,n] as described in the paper at
        http://www.diva-portal.org/smash/get/diva2:282089/fulltext01
    """
    return pow(2*alpha,0.75) * _sqrt(
            pow(4*alpha,l+m+n) / ( dfac(2*l-1) * dfac(2*m-1) * dfac(2*n-1) )
            )

def Sxyz(a,b,diffA,diffB,gamma):
    """
    Calculate the one dimensional overlap integral over two Gaussian functions
    divided by sqrt(pi)

    a,b: each an int
        The exponent of the polynom (e.g. x-directoim: (x-X)^a
    diffA, diffB: each a float
        The difference between the center of the combined Gaussian
        and the center of thje original Gaussian
    gamma: float
        The combined exponent
    """
    indices   = ((i,j) for i in xrange(0,a+1) for j in xrange(0,b+1) if (i+j)%2==0)
    result    = sum( (
        binomial(a,i) * binomial(b,j) * dfac(i+j-1) * pow(diffA,a-i) * pow(diffB,b-j) / pow(2*gamma,(i+j)*0.5)
        for i,j in indices
        ) )
    result /= _sqrt(gamma)
    return result
    

def S(A,B,alpha,beta,L1,L2):
    """
    Compute the overlap between two Cartesian Gaussian orbitals.
    
    A: list of 3 float
        Coordinates of center of first Gaussian
    B: list of 3 float
        Coordinates of center of second Gaussian
    alpha: float
        Exponential factor of first Gaussian
    beta: float
        Exponential factor of second Gaussian
    L1: list of 3 int
        The vector [l1,m1,n1] as described in the paper at
        http://www.diva-portal.org/smash/get/diva2:282089/fulltext01
    L2: list of 3 int
        Compare L1
    """
    gamma    = float(alpha+beta)
    eta      = alpha*beta/gamma
    P        = [(alpha*a + beta*b)/gamma for a,b in zip(A,B)]
    norm_2   = sum((a-b)*(a-b) for a,b in zip(A,B))
    EAB      = _exp(-eta*norm_2)
    iterator = ((a,b,Pi-Ai,Pi-Bi) for a,b,Ai,Bi,Pi in zip(L1,L2,A,B,P))
    if useFD:
        result   = reduce(float.__mul__, ( 
                                            SxyzPy(a,b,diffA,diffB,gamma) 
                                         for a,b,diffA,diffB in iterator 
                                         ) )
        result *= EAB * NormCoeffPy(alpha,*L1) * NormCoeffPy(beta,*L2)
    else:
        result   = reduce(float.__mul__, ( 
                                            Sxyz(a,b,diffA,diffB,gamma) 
                                         for a,b,diffA,diffB in iterator 
                                         ) )
        result *= EAB * normalization_coeff(alpha,*L1) * normalization_coeff(beta,*L2)
    return result

def Smatrix(basis,basis2=None):
    """
    Compute the overlap matrix of a given basis.

    basis: a list of [A,L,Prim]
           with
           A: a list of 3 floats
                The center of the contracted Cartesian Gaussian function
           L: a list of 3 ints
                The polynomial exponents of the contracted Cartesian Gaussian
           Prim: a list of [alpha,pre]
                with
                alpha: float
                    The exponential factor of the primitive Gaussian function
                pre: float
                    The contraction coefficient of the primitive Gaussian function
    basis2: same format as basis
        If not None, the overlap between basis and basis2 will be computed.
        Otherwise the self-overlap of basis will be computed.
    """
    if basis2 is None:
        basis2=basis
    #all matching brackets are alinged vertically
    #all for-statements are aligned with the closing bracket they belong to
    result = [    [   sum((
                           prefactorA*prefactorB * S(A,B,alpha,beta,L1,L2)
                       for alpha,prefactorA in PrimA for beta,prefactorB in PrimB
                       ))
                for B,L2,PrimB in basis2
                ] 
           for A,L1,PrimA in basis
           ]
    return result
    #this code does the same but is insanely slower:
    #matrix = [[0.0]*len(basis)]*len(basis)
    #i=0
    #for A,L1,PrimA in basis:
    #    j=0
    #    for B,L2,PrimB in basis:
    #        result = 0.0
    #        for alpha,prefactorA in PrimA:
    #            for beta,prefactorB in PrimB:
    #                result += prefactorA*prefactorB * S(A,B,alpha,beta,L1,L2)
    #        matrix[i][j] = result
    #        j+=1
    #    i+=1
    #return matrix

def normalize_basis(basis,Smat=None):
    if Smat is None:
        Smat = Smatrix(basis)
    correction = [1.0/_sqrt(Smat[i][i]) for i in xrange(len(Smat))]
    for i in xrange(len(basis)):
        for j in xrange(len(basis[i][2])):
            basis[i][2][j][1] *= correction[i]
    return [[Smat[i][j]*correction[i]*correction[j] for j in xrange(len(Smat))] for i in xrange(len(Smat))]

def overlap_lincomb(Smat,coefficients1,coefficients2=None):
    """
    Given the overlap matrix Smat and some coefficients, compute
    the overlap of two MOs expanded in terms of the basis that
    results in the given overlap matrix.

    If coefficients2 is not specified, the same coefficients are
    used twice.
    """
    if coefficients2 is None:
        coefficients2 = coefficients1
    if len(coefficients1) != len(Smat) or len(coefficients2) != len(Smat[0]):
        raise ValueError("Wrong dimensions for overlap compuation.")
    return sum((
                  ci*cj*Sij 
              for    ci,Si  in zip(coefficients1,Smat)    #outer loop
                 for cj,Sij in zip(coefficients2,Si)      #inner loop
              ))

def normalize_MOs(Smat,coefficients,occupations=None):
    """
    Normalize a molecular orbital defined by some coefficients to the
    given occupation. The overlap matrix for the current basis also needs
    to be specified.

    If coefficients is a nested list, the operation will be performed
    for every sub-list. If occupations is not None, it has to have an
    appropriate shape.
    """
    d = _depth(coefficients)
    if d==1:
        correction = _sqrt(1.0*occupations/overlap_lincomb(Smat,coefficients))
        for i in xrange(len(coefficients)):
            coefficients[i] *= correction
    elif d == 2:
        for j in xrange(len(coefficients)):
            correction = _sqrt(1.0*occupations[j]/overlap_lincomb(Smat,coefficients[j]))
            for i in xrange(len(coefficients[j])):
                coefficients[j][i] *= correction
    else:
        raise TypeError("You either have to specify a list of coefficients or a list of such lists. Depth of nested list is wrong.")
