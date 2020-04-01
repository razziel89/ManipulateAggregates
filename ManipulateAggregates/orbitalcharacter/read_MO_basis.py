"""Treat Cartesian Gaussian basis set information.

High-level module that allows for reading in basis set information and
transforming it in a format suitable for further processing.
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
import copy
import math
import sys

import logging

logger = logging.getLogger(__name__)
try:
    from ..collection.read import read_molden
except ImportError:
    logger.warning("Could not import ..collection.read.read_molden")


def _get_MOs_occupation(MOs, occ_func):
    """Filter molecular orbitals with respect to their occupations and separate
       them into alpha and beta spins. Also determine whether the alpha HOMO is
       higher in energy than the beta HOMO. Also determine occupations.
    """
    MOsalpha = [mo[3] for mo in MOs if mo[1].lower() == "alpha"]
    MOsbeta = [mo[3] for mo in MOs if mo[1].lower() == "beta"]
    OCCsalpha = [mo[2] for mo in MOs if mo[1].lower() == "alpha"]
    OCCsbeta = [mo[2] for mo in MOs if mo[1].lower() == "beta"]
    if len(MOsalpha) > 0 and len(MOsbeta) > 0:
        Esalpha = [mo[0] for mo in MOs if mo[1].lower() == "alpha" and occ_func(mo[2])]
        Esbeta = [mo[0] for mo in MOs if mo[1].lower() == "beta" and occ_func(mo[2])]
        max_energy_alpha = max(Esalpha)
        max_energy_beta = max(Esbeta)
        IdxHOMOalpha = len(Esalpha) - 1
        IdxHOMObeta = len(Esbeta) - 1
        high_alpha = max_energy_alpha >= max_energy_beta
    else:
        Es = [mo[0] for mo in MOs if occ_func(mo[2])]
        IdxHOMOalpha = len(Es) - 1
        IdxHOMObeta = IdxHOMOalpha
        high_alpha = True
    if len(MOsalpha) > 0 or len(MOsbeta) > 0:
        return (
            high_alpha,
            (MOsalpha, OCCsalpha),
            (MOsbeta, OCCsbeta),
            (IdxHOMOalpha, IdxHOMObeta),
        )
    else:
        raise ValueError(
            "No molecular orbitals fulfilling the occupation function defined for either spin."
        )


def _gen_basis_from_GTO(Atoms, GTO, sorting="none"):
    """Get a basis of contracted Cartesian Gaussian functions from contracted shells.

    The shells are expanded in terms of their angular momenta.

    Args:
        Atoms: (list of [ignored, [3 floats]]) "ignored" is ignored. The three
            floats are the Cartesian coordinates of the atoms.
        GTO: (list of [shelltype,prefactor,nr_prim,prim]) "shelltype" is a letter
            declaring the type of shell. "prefactor" is ignored by the molden
            standard. "nr_prim" is an integer number declaring the number of
            primitives in the shell. "prim" is a list of lists of 2 floats
            declaring the exponential and prefactor of the primitives of this
            shell.
        sorting: can be 'none' (s,p,d,f orbitals for each atom one after the other)
            or 'spdf' (have an ss block, then an sp-block, etc.)

    Raises:
        ValueError.

    Yiels:
        Contracted Cartesian Gaussian orbitals in the form of
        (x,y,z),(n,l,m),prim where x,y and z are the Cartesian coordinates of
        the center of the atomic orbital. n,l and m are the angular momentum
        coefficients (integers) and prim is what's mentioned for GTO.
    """
    GTO = [gto[1] for gto in GTO]
    if sorting.lower() == "spdf":
        for wanttype in ["s", "p", "d"]:
            for (element, (x, y, z)), gto in zip(Atoms, GTO):
                # prefactor is deliberately being ignored because it is no longer
                # used by the molden file standard
                for shelltype, prefactor, nr_prim, prim in gto:
                    if shelltype == wanttype:
                        if shelltype == "s":
                            yield (x, y, z), (0, 0, 0), prim
                        elif shelltype == "p":
                            yield (x, y, z), (1, 0, 0), prim
                            yield (x, y, z), (0, 1, 0), prim
                            yield (x, y, z), (0, 0, 1), prim
                        elif shelltype == "d":
                            yield (x, y, z), (2, 0, 0), prim
                            yield (x, y, z), (0, 2, 0), prim
                            yield (x, y, z), (0, 0, 2), prim
                            yield (x, y, z), (1, 1, 0), prim
                            yield (x, y, z), (1, 0, 1), prim
                            yield (x, y, z), (0, 1, 1), prim
                        elif shelltype == "f":
                            yield (x, y, z), (3, 0, 0), prim
                            yield (x, y, z), (0, 3, 0), prim
                            yield (x, y, z), (0, 0, 3), prim
                            yield (x, y, z), (2, 1, 0), prim
                            yield (x, y, z), (2, 0, 1), prim
                            yield (x, y, z), (1, 2, 0), prim
                            yield (x, y, z), (0, 2, 1), prim
                            yield (x, y, z), (1, 0, 2), prim
                            yield (x, y, z), (0, 1, 2), prim
                            yield (x, y, z), (1, 1, 1), prim
                        else:
                            raise ValueError(
                                "Only s,p,d and f-type shells are supported."
                            )
    elif sorting.lower() == "none":
        for (element, (x, y, z)), gto in zip(Atoms, GTO):
            # prefactor is deliberately being ignored because it is no longer
            # used by the molden file standard
            for shelltype, prefactor, nr_prim, prim in gto:
                if shelltype == "s":
                    yield (x, y, z), (0, 0, 0), prim
                elif shelltype == "p":
                    yield (x, y, z), (1, 0, 0), prim
                    yield (x, y, z), (0, 1, 0), prim
                    yield (x, y, z), (0, 0, 1), prim
                elif shelltype == "d":
                    yield (x, y, z), (2, 0, 0), prim
                    yield (x, y, z), (0, 2, 0), prim
                    yield (x, y, z), (0, 0, 2), prim
                    yield (x, y, z), (1, 1, 0), prim
                    yield (x, y, z), (1, 0, 1), prim
                    yield (x, y, z), (0, 1, 1), prim
                elif shelltype == "f":
                    yield (x, y, z), (3, 0, 0), prim
                    yield (x, y, z), (0, 3, 0), prim
                    yield (x, y, z), (0, 0, 3), prim
                    yield (x, y, z), (1, 2, 0), prim
                    yield (x, y, z), (2, 1, 0), prim
                    yield (x, y, z), (2, 0, 1), prim
                    yield (x, y, z), (1, 0, 2), prim
                    yield (x, y, z), (0, 1, 2), prim
                    yield (x, y, z), (0, 2, 1), prim
                    yield (x, y, z), (1, 1, 1), prim
                else:
                    raise ValueError("Only s,p,d and f-type shells are supported.")
    else:
        raise ValueError("Sorting must be None or SPD")


# default function to check for occupations (consider everything as occupied)
global OCC_FUNC
OCC_FUNC = lambda o: True
# all variable names that contain "alpha" or "beta" are spin-polarized values
def get_MOs_and_basis(
    filename,
    occ_func=OCC_FUNC,
    filetype="molden",
    spins="alpha",
    copy_non_present=True,
    msave=None,
    sorting="none",
    alpha_high_energy=False,
):
    """Read information about molecular orbitals and the basis from a file.

    Basis functions are contracted Cartesian Gaussian functions.

    Args:
        filename: (string) the name of the file from which to read in the data
        occ_func: (function, returns boolean) a function that selects orbitals
            based on their occupation numbers. The indices that are returned
            correspond to the HOMOs of the desired spins. Give "lambda o: o>0.5"
            (without the quotes) if you want to consider only orbitals that are
            more occupied than 0.5
        filetype: (string) the type of file from which to read the data
            Possible options are "molden"
        spins: (str) declare which spins to return. Possible values are "alpha",
            "beta" and "both" depending on what is desired
        copy_non_present: (bool) if a requested spin is not present in the
            file, replace it by the present one. That way declaring 'both' for
            spins will return two lists even if no spin 'beta' or 'alpha'
            is given.
        msave: (dictionary) if a dictionary is given, add data from the
            internal variable m (return value of
            ManipulateAggregates.collection.read.read_molden. This makes
            printing out another molden file using the same basis
        sorting: (string) "none" (s,p,d,f orbitals for each atom one after the
            other) or "spdf" (have an ss block, then an sp-block, etc. in
            the basis)
        alpha_high_energy: (bool) if True, make it so that, if the selected
            orbitals have different energies, depending on their spins, the
            orbital of alpha spin is that one with the higher energy.
    """

    if spins == "both":
        want_alpha = True
        want_beta = True
    elif spins == "alpha":
        want_alpha = True
        want_beta = False
    elif spins == "beta":
        want_alpha = False
        want_beta = True
    else:
        raise ValueError("Unknown spin requested.")

    if filetype == "molden":
        m = read_molden(
            filename,
            positions=True,
            elementnames=True,
            GTO=True,
            GTO_coefficients=True,
            GTO_nr_primitives=False,
            MO=True,
            MO_coefficients=True,
        )
    else:
        raise ValueError("Unsupported filetype specified.")

    (
        high_alpha,
        (MOsalpha, OCCsalpha),
        (MOsbeta, OCCsbeta),
        (IdxHOMOalpha, IdxHOMObeta),
    ) = _get_MOs_occupation(m["MO"], occ_func)

    # If the user wants the alpha spin to be the energetically higher one (with respect to
    # the orbitals determined by occ_func) but alpha is not, swap both spins. This implies
    # that both alpha and beta spins are present, i.e., len(MOsalpha)>0 and len(MOsbeta)>0
    # so that the "want_alpha/beta" part will not be triggered if a swap occurred
    if alpha_high_energy and not (high_alpha):
        (MOsbeta, OCCsbeta, IdxHOMObeta), (MOsalpha, OCCsalpha, IdxHOMOalpha) = (
            (MOsalpha, OCCsalpha, IdxHOMOalpha),
            (MOsbeta, OCCsbeta, IdxHOMObeta),
        )

    basis = list(_gen_basis_from_GTO(m["positions"], m["GTO"]))

    # one of the two spins is defintely present, otherwise _get_MOs_occupation raised an error
    if want_alpha and len(MOsalpha) == 0:
        if copy_non_present:
            MOsalpha = copy.deepcopy(MOsbeta)
            # if one spin is duplicated, the occupations have to be halved
            OCCsbeta = [0.5 * occ for occ in OCCsbeta]
            OCCsalpha = copy.deepcopy(OCCsbeta)
    if want_beta and len(MOsbeta) == 0:
        if copy_non_present:
            MOsbeta = copy.deepcopy(MOsalpha)
            # if one spin is duplicated, the occupations have to be halved
            OCCsalpha = [0.5 * occ for occ in OCCsalpha]
            OCCsbeta = copy.deepcopy(OCCsalpha)

    if not want_alpha:
        del MOsalpha
        MOsalpha = None
        del OCCsalpha
        OCCsalpha = None
        del IdxHOMOalpha
        IdxHOMOalpha = None
    if not want_beta:
        del MOsbeta
        MOsbeta = None
        del OCCsbeta
        OCCsbeta = None
        del IdxHOMObeta
        IdxHOMObeta = None

    if msave is not None:
        for key in m:
            msave[key] = copy.deepcopy(m[key])

    return (
        basis,
        (MOsalpha, OCCsalpha),
        (MOsbeta, OCCsbeta),
        (IdxHOMOalpha, IdxHOMObeta),
    )
