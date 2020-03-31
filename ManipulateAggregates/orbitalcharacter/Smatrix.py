"""Compute overlaps of basis functions and molecular orbitals.

Taken from http://www.mathematica-journal.com/2012/02/evaluation-of-gaussian-molecular-integrals/
and altered a bit before transforming to Python code. All thanks goes to the authors.
"""

# This file is part of ManipulateAggregates.
#
# Copyright (C) 2016 by Torsten Sachse
#
# ManipulateAggregates is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ManipulateAggregates is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from collections import Sequence
from itertools import chain, count
import math

_exp = math.exp
_sqrt = math.sqrt

try:
    # these are faster C++ routines that do the same thing
    from FireDeamon import normalization_coefficient as NormCoeffPy
    from FireDeamon import Sxyz as SxyzPy

    useFD = True
except ImportError:
    useFD = False
if not useFD:
    try:
        from scipy.special import binom as binomial
    except ImportError:
        from math import factorial

        binomial = lambda n, k: factorial(n) / (factorial(n - k) * factorial(k))


def _depth(seq):
    """
    Calculate the depth of a list (seq). Taken from
    http://stackoverflow.com/questions/6039103/counting-deepness-or-the-deepest-level-a-nested-list-goes-to
    """
    seq = iter(seq)
    try:
        for level in count():
            seq = chain([next(seq)], seq)
            seq = chain.from_iterable(s for s in seq if isinstance(s, Sequence))
    except StopIteration:
        return level


def _dfac(n):
    """
    Calculate the double factorial (not found in module "math")

    n: int
        Integer number of which to calculate the double factorial
    """
    x = n if n >= 0 else 1
    for y in range(n - 2, 1, -2):
        x *= y
    return x


def normalization_coeff(alpha, l, m, n):
    """Calculate the normalization coefficient of a 3D Cartesian Gaussian function times pi^0.75.

    The list [l,m,n] as described in the paper at http://www.diva-portal.org/smash/get/diva2:282089/fulltext01

    Returns:
        the normalization coefficient.

    Args:
        alpha: (float) exponential factor of the Gaussian
        l: (int) angular factor
        m: (int) angular factor
        n: (int) angular factor
    """
    return pow(2 * alpha, 0.75) * _sqrt(
        pow(4 * alpha, l + m + n)
        / (_dfac(2 * l - 1) * _dfac(2 * m - 1) * _dfac(2 * n - 1))
    )


def Sxyz(a, b, diffA, diffB, gamma):
    """Calculate the one dimensional overlap integral over two Gaussian functions divided by sqrt(pi).

    Returns:
        overlap of two Gaussian functions

    Args:
        a: (int) exponent of the first Cartesian prefactor (x-x_a0)^a
        b: (int) exponent of the second Cartesian prefactor (x-x_b0)^b
        diffA: (float) the difference between the center of the combined
            Gaussian and the center of the first original Gaussian
        diffB: (float) the difference between the center of the combined
            Gaussian and the center of the second original Gaussian
        gamma: (float) the combined exponent
    """
    indices = (
        (i, j) for i in range(0, a + 1) for j in range(0, b + 1) if (i + j) % 2 == 0
    )
    result = sum(
        (
            binomial(a, i)
            * binomial(b, j)
            * _dfac(i + j - 1)
            * pow(diffA, a - i)
            * pow(diffB, b - j)
            / pow(2 * gamma, (i + j) * 0.5)
            for i, j in indices
        )
    )
    result /= _sqrt(gamma)
    return result


def S(A, B, alpha, beta, L1, L2):
    """Compute the overlap between two Cartesian Gaussian orbitals.
    
    Args:
        A: (list of 3 floats) coordinates of center of first Gaussian
        B: (list of 3 floats) coordinates of center of second Gaussian
        alpha: (float) exponential factor of first Gaussian
        beta: (float) exponential factor of second Gaussian
        L1: (list of 3 int) The vector [l1,m1,n1] as described in the paper
            at http://www.diva-portal.org/smash/get/diva2:282089/fulltext01
        L2: (list of 3 int) see L1

    Returns:
        overlap between two Cartesian Gaussian orbitals.
    """
    gamma = float(alpha + beta)
    eta = alpha * beta / gamma
    P = [(alpha * a + beta * b) / gamma for a, b in zip(A, B)]
    norm_2 = sum((a - b) * (a - b) for a, b in zip(A, B))
    EAB = _exp(-eta * norm_2)
    iterator = ((a, b, Pi - Ai, Pi - Bi) for a, b, Ai, Bi, Pi in zip(L1, L2, A, B, P))
    if useFD:
        hereSxyz = SxyzPy
        hereNormCoeff = NormCoeffPy
    else:
        hereSxyz = Sxyz
        hereNormCoeff = normalization_coeff
    result = 1.0
    for a, b, diffA, diffB in iterator:
        result *= hereSxyz(a, b, diffA, diffB, gamma)
    result *= EAB * hereNormCoeff(alpha, *L1) * hereNormCoeff(beta, *L2)
    return result


def Smatrix(basis, basis2=None):
    """Compute the overlap matrix of a given basis.

    Args:
        basis: (a list of [A,L,Prim])
            object of ManipulateAggregates.orbitalcharacter.density_on_grid.basis
        basis2: (same format as basis) if not None, the overlap between basis
            and basis2 will be computed. Otherwise the self-overlap of basis will
            be computed.

    Returns:
        a square matrix that is the overlap between the shells in the basis
    """
    if basis2 is None:
        basis2 = basis
    # all matching brackets are alinged vertically
    # all for-statements are aligned with the closing bracket they belong to
    result = [
        [
            sum(
                (
                    prefactorA * prefactorB * S(A, B, alpha, beta, L1, L2)
                    for alpha, prefactorA in PrimA
                    for beta, prefactorB in PrimB
                )
            )
            for B, L2, PrimB in basis2
        ]
        for A, L1, PrimA in basis
    ]
    return result
    # this code does the same but is insanely slower:
    # matrix = [[0.0]*len(basis)]*len(basis)
    # i=0
    # for A,L1,PrimA in basis:
    #    j=0
    #    for B,L2,PrimB in basis:
    #        result = 0.0
    #        for alpha,prefactorA in PrimA:
    #            for beta,prefactorB in PrimB:
    #                result += prefactorA*prefactorB * S(A,B,alpha,beta,L1,L2)
    #        matrix[i][j] = result
    #        j+=1
    #    i+=1
    # return matrix


def normalize_basis(basis, Smat=None):
    """Normalize a basis of Cartesian gaussian orbitals.

    Args:
        basis: (a list of [A,L,Prim])
               object of ManipulateAggregates.orbitalcharacter.density_on_grid.basis
        Smat: (a square matrix of floats) the overlap between the shells in the
            basis. If it is None, a suitable overlap matrix will be computed.

    Returns:
        The overlap matrix for the normalized basis. The given basis has
        already been adjusted.
    """
    if Smat is None:
        Smat = Smatrix(basis)
    correction = [1.0 / _sqrt(Smat[i][i]) for i in range(len(Smat))]
    for i in range(len(basis)):
        for j in range(len(basis[i][2])):
            basis[i][2][j][1] *= correction[i]
    return [
        [Smat[i][j] * correction[i] * correction[j] for j in range(len(Smat))]
        for i in range(len(Smat))
    ]


def overlap_lincomb(Smat, coefficients1, coefficients2=None):
    """Compute the overlap of molecular orbitals.

    Given the overlap matrix Smat and some coefficients, compute the overlap of
    two MOs expanded in terms of the basis that results in the given overlap
    matrix.

    Args:
        Smat: (a square matrix of floats) the overlap between the shells in the
            basis. If it is none, a suitable overlap matrix will be computed.
        coefficients1: (lists of floats) these coefficients define a molecular
            orbital in the basis whose overlap matrix is also given.
        coefficients2: (list of floats) if provided, the overlap between the MO
            specified by coefficients1 and this parameter is computed.

    Returns:
        the overlap of the molecular orbitals
    """
    if coefficients2 is None:
        coefficients2 = coefficients1
    if len(coefficients1) != len(Smat) or len(coefficients2) != len(Smat[0]):
        raise ValueError("Wrong dimensions for overlap compuation.")
    return sum(
        (
            ci * cj * Sij
            for ci, Si in zip(coefficients1, Smat)  # outer loop
            for cj, Sij in zip(coefficients2, Si)  # inner loop
        )
    )


def normalize_MOs(Smat, coefficients, occupations=None):
    """Normalize a molecular orbital defined by some coefficients to the given occupation.

    Args:
        Smat: (a square latrix of floats) the overlap between the shells in the
            basis. If it is none, a suitable overlap matrix will be computed.
        coefficients: (list of floats or a list of lists of floats) these
            coefficients define molecular orbital in the basis. If only one
            orbital is given, it is normalized.
        occupations: (list of floats) if provided the molecular orbitals will
            be normalized to the occupations in this list. If only one orbital
            is given, it is normalized to the number given here. everything is
            normalized to 1 by default.

    Raises:
        TypeError.
    """
    d = _depth(coefficients)
    if d == 1:
        if occupations is None:
            occupations = 1.0
        correction = _sqrt(1.0 * occupations / overlap_lincomb(Smat, coefficients))
        for i in range(len(coefficients)):
            coefficients[i] *= correction
    elif d == 2:
        if occupations is None:
            occupations = [1.0] * len(coefficients)
        for j in range(len(coefficients)):
            correction = _sqrt(
                1.0 * occupations[j] / overlap_lincomb(Smat, coefficients[j])
            )
            for i in range(len(coefficients[j])):
                coefficients[j][i] *= correction
    else:
        raise TypeError(
            "You either have to specify a list of coefficients or a list of such lists. Depth of nested list is wrong."
        )
