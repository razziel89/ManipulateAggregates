"""This submodule agregates functions that can be applied to QM orbitals.
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
import sys
import copy

import logging

logger = logging.getLogger(__name__)
try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")

try:
    from . import read_MO_basis
except ImportError:
    logger.warning("Could not import .read_MO_basis")
try:
    from . import Smatrix
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
    logger.warning(
        "Could not import ..orbitalcharacter.read_MO_basis.get_MOs_and_basis"
    )
try:
    from ..orbitalcharacter.Smatrix import (
        Smatrix,
        normalize_basis,
        overlap_lincomb,
        normalize_MOs,
    )
except ImportError:
    logger.warning(
        "Could not import Smatrix, normalize_basis, overlap_lincomb or normalize_MOs from ..orbitalcharacter.Smatrix"
    )

try:
    from FireDeamon import ElectronDensityPy, InitializeGridCalculationOrbitalsPy
except ImportError:
    logger.warning(
        "Could not import ElectronDensityPy or InitializeGridCalculationOrbitalsPy from FireDeamon"
    )
try:
    from FireDeamon import (
        InitializeGridCalculationOrbitalsPy,
        ElectrostaticPotentialOrbitalsPy,
        ElectrostaticPotentialPy,
    )
except ImportError:
    logger.warning(
        "Could not import InitializeGridCalculationOrbitalsPy, ElectrostaticPotentialOrbitalsPy or ElectrostaticPotentialPy from FireDeamon"
    )

# all orbitals with a higher occupation than this are considered to be occupied
MINOCC = 0.0000001
# convertion factor from Bohr (atomic units) to Angstroms
BOHRTOANG = 0.529177249


def expand_total_wavefunction(MOs, Smat, normalize=False):
    """Compute the total wavefunction, i.e., the sum over all molecular orbitals.

    Args:
        MOs: (list of lists of floats) coefficients describing the molecular orbitals
        Smat: (square matrix of floats) matrix describing the overlap
            between the shels of the basis used to obtain the molecular orbitals
        normalize: (bool) whether or not to make it so that the overlap of the
            returned total wavefunction with itself shall be 1

    Returns:
        a list of floats describing the total wavefunction in the given basis
    """
    Bsize = len(Smat)
    MOsize = len(MOs)
    RHS = [
        sum((coeff * Skj for mo in MOs for Skj, coeff in zip(Sk, mo))) for Sk in Smat
    ]
    # I want to have: |P> = SUM (d_k*|phi_k>) for i in inteval [0,N]
    # with: P: total wavefunction, d_k: linear combination coefficient
    #      S_ik: element of the overlap matrix i.e. <phi_i|phi_k> (symmetric)
    #      N: nr. basis functions
    #      |phi_k>: k-th basis function
    # Hence, <phi_i|P> = SUM (d_k <phi_i|phi_k>) for k in interval [0,N]
    #                 = SUM (d_k S_ik) for k in interval [0,N]
    # With v = [<phi_1|P>,<phi_2|P>,...,<phi_N|P>] and d = [d_1,d_2,...,d_N] (Python lists as vectors)
    #     follows v = S dot d (dot: matrix product)
    #     and hence the equation S*d=v has to be solved for v
    # solves S*d=RHS where d is the vector containing the coefficients of interest
    result = numpy.linalg.solve(Smat, RHS)
    if normalize:
        result /= sqrt(overlap_lincomb(Smat, result))
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
    print(
        "DEBUG: reading in basis from molden file and applying corrections for limited precision read",
        file=sys.stderr,
    )
    # read in the molden file and extract spin-polarized MO information from it
    # also read in basis information
    occ_func = lambda o: o > MINOCC
    (
        basis,
        (allMOsalpha, allOCCsalpha),
        (allMOsbeta, allOCCsbeta),
        (IdxHOMOalpha, IdxHOMObeta),
    ) = get_MOs_and_basis(
        filename,
        filetype="molden",
        spins="both",
        alpha_high_energy=True,
        occ_func=occ_func,
    )
    # determine occupied orbitals
    MOsalpha = allMOsalpha[: IdxHOMOalpha + 1]
    MOsbeta = allMOsbeta[: IdxHOMObeta + 1]
    OCCsalpha = allOCCsalpha[: IdxHOMOalpha + 1]
    OCCsbeta = allOCCsbeta[: IdxHOMObeta + 1]
    Smat = Smatrix(basis)
    Smat = normalize_basis(basis, Smat)  # after this, basis will be normalized
    copy_beta = MOsalpha == MOsbeta
    normalize_MOs(Smat, MOsalpha, occupations=OCCsalpha)
    if copy_beta:
        MOsbeta = copy.deepcopy(MOsalpha)
    else:
        normalize_MOs(Smat, MOsbeta, occupations=OCCsbeta)
    return basis, Smat, (MOsalpha, MOsbeta), (OCCsalpha, OCCsbeta)


def _correlation(array1, array2):
    """Return the correlation between two NumPy arrays."""
    temparray1 = array1 - numpy.mean(array1)
    temparray2 = array2 - numpy.mean(array2)
    dot1 = numpy.dot(temparray1, temparray1)
    dot2 = numpy.dot(temparray2, temparray2)
    temparray1 /= numpy.sqrt(dot1)
    temparray2 /= numpy.sqrt(dot2)
    return numpy.dot(temparray1, temparray2)


def _similarity(diffdens, MOdens, type=0, name=False):
    """Compute the similarity between a difference density
    and the density of a molecular orbital.
    """
    if type == 0:
        # This should be as close to 1 as possible
        # but I guess the other two are a better
        # meassure
        result = _correlation(diffdens, MOdens)
        rname = "Correlation"
    elif type == 1:
        # This should be as close to 1 as possible
        # which means that a lot of density is being
        # taken from the MO
        temparray1 = numpy.copy(diffdens)
        temparray2 = MOdens
        temparray1[temparray1 < 0.0] = 0.0
        result = _correlation(temparray1, temparray2)
        rname = "Correlation Positive"
    elif type == 2:
        # This should be as close to 0 as possible
        # which means that no density is being transferred
        # into the MO
        temparray1 = numpy.copy(diffdens)
        temparray2 = MOdens
        temparray1[temparray1 > 0.0] = 0.0
        temparray1 *= -1
        result = _correlation(temparray1, temparray2)
        rname = "Correlation Negative"
    else:
        raise ValueError("Wrong type of similarity measure given.")
    if name:
        return result, rname
    else:
        return result


def density_overlap(density_1, density_2):
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
    print("DEBUG: started computation of density correlation", file=sys.stderr)
    # read in all dx files
    data1 = numpy.array(
        read_dx(density_1, density=True, silent=True, grid=False, gzipped=True)["data"]
    )
    data2 = numpy.array(
        read_dx(density_2, density=True, silent=True, grid=False, gzipped=True)["data"]
    )
    print("DEBUG: reading dx-files done", file=sys.stderr)
    if data1.shape != data2.shape:
        raise ValueError(
            "Both dx files contain grids with a different number of points."
        )
    corrvalue = _similarity(data1, data2, 0, name=False)
    print("DEBUG: computed overlap between densities", file=sys.stderr)
    print("Overlap: %8.4e" % (corrvalue))
    print("DEBUG: done computation of density correlation", file=sys.stderr)


def difference_density(density_1, density_2, dxdiffdens, compress=False, factor=1.0):
    """Just compute the difference between two densities given by their dx-files.

    Args:
        density_1: (string) name of the dx-file that contains the first density
        density_2: (string) name of the dx-file that contains the second density
        dxdiffdens: (string) name of the dx-file that shall contain the difference density.
        compress: (bool) whether or not the difference density shall be written
            in gzipped format or not.
        factor: (float) this factor is multiplied with the second density
            before computing the difference
    """
    print("DEBUG: started computation of difference density", file=sys.stderr)
    apply_func_density(
        density_1,
        density_2,
        dxdiffdens,
        compress=compress,
        func=lambda d1, d2: d1 - factor * d2,
        verbose=True,
    )
    print("DEBUG: done computation of difference density", file=sys.stderr)


global DEFAULT_FUNC
# default function for ManipulateAggregates.orbitalcharacter.apply_func_density
DEFAULT_FUNC = lambda d1, d2: d1


def apply_func_density(
    density_1, density_2, outdens, compress=False, func=DEFAULT_FUNC, verbose=False
):
    """Apply a given function to two densities and write result to file.

    If density_2 is None, apply the given function only to the first
    density. If no function is provided, only the first density is written out.

    Bug:
        If both dx-files define unequal grids that have the same number of
        points, the correlation is still computed but the results are not the
        actual correlations.

    Raises:
        ValueError

    Args:
        density_1: (string) - name of the dx-file that contains the first density
        density_2: (string) - name of the dx-file that contains the second density
        outdens: (string) - name of the dx-file that shall contain the output density
        compress: (bool) - whether or not the difference density shall be written in
            gzipped format or not
        func: (function of 2 variables applicable to numpy arrays) - How to obtain the
            new density when given the old ones
        verbose: (bool) - give progress updates

    """
    header = {}
    # read in all dx files
    data1 = numpy.array(
        read_dx(
            density_1,
            density=True,
            silent=True,
            grid=False,
            header_dict=header,
            gzipped=True,
        )["data"]
    )
    if density_2 is not None:
        data2 = numpy.array(
            read_dx(density_2, density=True, silent=True, grid=False, gzipped=True)[
                "data"
            ]
        )
    if verbose:
        print("DEBUG: reading dx-files done", file=sys.stderr)
    if density_2 is not None:
        if data1.shape != data2.shape:
            raise ValueError(
                "Both dx files contain grids with a different number of points."
            )
        newdens = func(data1, data2)
    else:
        newdens = func(data1)
    if verbose:
        print(
            "DEBUG: computed new density, sum: %8.4f, sum over abs: %8.4f"
            % (numpy.sum(diffdens), numpy.sum(numpy.fabs(diffdens))),
            file=sys.stderr,
        )
    print_dx_file(
        outdens,
        header["counts_xyz"],
        header["org_xyz"],
        header["delta_x"],
        header["delta_y"],
        header["delta_z"],
        newdens,
        gzipped=compress,
    )
    if verbose:
        print("DEBUG: wrote new density", file=sys.stderr)


def _postprocess_multiple(
    total_1, total_2, MOalpha_1, MObeta_1, MOalpha_2, MObeta_2, dir="", type="kation"
):
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
    print("DEBUG: started postprocessing", file=sys.stderr)
    if len(dir) > 0:
        if not dir.endswith("/"):
            dir += "/"
    # read in all dx files
    # they have been normalized to the number of electrons
    header = {}
    data1 = numpy.array(
        read_dx(
            total_1,
            density=True,
            silent=True,
            grid=False,
            header_dict=header,
            gzipped=True,
        )["data"]
    )
    data2 = numpy.array(
        read_dx(total_2, density=True, silent=True, grid=False, gzipped=True)["data"]
    )
    print("DEBUG: reading dx-files done", file=sys.stderr)
    nr_electrons_1 = int(round(numpy.sum(data1)))
    nr_electrons_2 = int(round(numpy.sum(data2)))
    if type == "kation":
        kation = True
        diff_to_neut = +1
        prefix = "kat"
    elif type == "anion":
        kation = False
        diff_to_neut = -1
        prefix = "an"
    else:
        raise ValueError("Wrong type of molecule comparison given.")
    if nr_electrons_1 == nr_electrons_2 + diff_to_neut:
        neut_total = data1
        ion_total = data2
        if type == "kation":
            MOalpha = numpy.array(
                read_dx(MOalpha_1, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
            MObeta = numpy.array(
                read_dx(MObeta_1, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
        else:
            MOalpha = numpy.array(
                read_dx(MOalpha_2, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
            MObeta = numpy.array(
                read_dx(MObeta_2, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
        nr_electrons_neut = nr_electrons_1
    elif nr_electrons_1 == nr_electrons_2 - diff_to_neut:
        neut_total = data2
        ion_total = data1
        if type == "kation":
            MOalpha = numpy.array(
                read_dx(MOalpha_2, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
            MObeta = numpy.array(
                read_dx(MObeta_2, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
        else:
            MOalpha = numpy.array(
                read_dx(MOalpha_1, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
            MObeta = numpy.array(
                read_dx(MObeta_1, density=True, silent=True, grid=False, gzipped=True)[
                    "data"
                ]
            )
        nr_electrons_neut = nr_electrons_2
    else:
        raise ValueError(
            "Both dx files contain data about molecules that do not differ in exactly one electron."
        )
    print("DEBUG: determined %sion and neutral molecule" % (prefix), file=sys.stderr)
    if neut_total.shape != ion_total.shape:
        raise ValueError(
            "Both dx files contain grids with a different number of points."
        )
    diffdens = (neut_total - ion_total) * diff_to_neut
    print(
        "DEBUG: computed difference density, sum: %8.4f, sum over abs: %8.4f"
        % (numpy.sum(diffdens), numpy.sum(numpy.fabs(diffdens))),
        file=sys.stderr,
    )
    print_dx_file(
        dir + "diff_%sion.dx" % (prefix),
        header["counts_xyz"],
        header["org_xyz"],
        header["delta_x"],
        header["delta_y"],
        header["delta_z"],
        diffdens,
        gzipped=True,
    )
    print("DEBUG: wrote dx-file for difference density", file=sys.stderr)
    otypes = 3
    overlap_alpha = tuple(
        _similarity(diffdens, MOalpha, t, name=True) for t in range(otypes)
    )
    overlap_beta = tuple(
        _similarity(diffdens, MObeta, t, name=True) for t in range(otypes)
    )
    print(
        "DEBUG: computed overlap between difference density and HOMO densities (both spins)",
        file=sys.stderr,
    )
    for (a, na), (b, nb) in zip(overlap_alpha, overlap_beta):
        print("Type %15s: Overlap alpha/beta: %8.4e /%8.4e" % (na, a, b))
    print("DEBUG: done postprocessing", file=sys.stderr)
