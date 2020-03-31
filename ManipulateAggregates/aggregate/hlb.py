# -*- coding: utf-8 -*-
"""Calculate the HLB-value using multiple methods.

There are two methods supported by this subsubmodule:
  1. simple: Scan the molecule with a plane that is allowed to cut midway
    through bonds with a maximum angle of 30° to the bond. Assign positive values
    to all hetero atoms and negative values to C and H. Try to maximize the
    number of heteroatoms on one side and the number of C and H on the other
    side.
  2. complex: Scan the molecule with a plane that is allowed to cut midway
    through bonds with a maximum angle of 30° to the bond. Try to maximize the
    amount of positive electrostatic potential on one side and try to minimze
    the amount of negative electrostatic potential on the other side of the
    plane.

I highly recommend using the simple method as gave values closer to measured
ones, up to now.
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
import logging

logger = logging.getLogger(__name__)

from math import pi, cos, sin, copysign, asin, acos
import sys
import re
import time

try:
    from numpy import (
        array,
        linalg,
        dot,
        arange,
        array_equal,
        inf,
        cross,
        sqrt,
        invert,
        average,
        amax,
        amin,
        fabs,
        concatenate,
    )
    from numpy import max
    from numpy.linalg import norm
except ImportError as e:
    logger.warning(
        "Could not import from numpy: array, linalg, dot, arange, array_equal, inf, cross, "
        + "sqrt, invert, average, amax, amin, fabs, concatenate, max, linalg.norm"
    )

global MASSES_DICT, WEIGHTS_DICT
# this dictionary stores all masses of the atoms to be used in relative units
MASSES_DICT = {
    "H": 1.0008,
    "O": 16.0,
    "S": 32.06,
    "N": 14.01,
    "F": 19.0,
    "Cl": 35.45,
    "C": 12.01,
    "Ni": 58.69,
}
# this dictionary stores all weights to judge which atom is rather hydrophilic (high
# positive value) and which is not
WEIGHTS_DICT = {"H": 0, "O": 1, "S": 1, "N": 1, "F": 1, "Cl": 1, "C": -0.4, "Ni": 1}


def part_molecule(structure, normal_vector, coordinate):
    """Part a molecule by a plane (in Hessian normal form).

    Returns:
        a list of indices of all elements of structure that are on the side
        of the plane (defined by normal_vector and coordinate) to which
        the normal vector points

    Args:
        structure: (list of vectors) Cartesian coordinates that shall be separated by
            the plane
        normal_vector: (list of 3 floats) normal vector for Hessian normal form
        coordinate: (list of 3 floats) point in the plane for Hessian normal form
    """
    return [
        i
        for i in range(0, len(structure))
        if dot(array(structure[i]) - array(coordinate), array(normal_vector)) >= 0
    ]


def list_min(l, key=lambda x: x):
    """Get the minimum element of a list.

    Returns:
        the minimum element of a list

    Args:
        l: (list) list whose minimum element shall be determined
        key: (function of 1 argument) this function is applied to every element
            of l to get an element that the "min" function understands
    """
    li = [key(le) for le in l]
    return li.index(min(li))


def pseudo_energy_simple(structure, names, normal_vector, point_in_plane):
    """Use the simple method to determine the optimal location of the air-water interface.

    With the given structure, calculate the "energy" associated with it if the
    plane, given in its Hessian normal form, cuts through the molecule. The
    value is high and positive, if atoms with a positive weight are on the side
    where the normal vector points and atoms with negative weights are on the
    opposite side of the plane.

    Returns:
        the pseudo energy associated with the current system (i.e., the plane
        defined by normal_vector and point_in_plane cutting the molecule
        defined by structure and names)

    Args:
        structure: (a list of 3-element vectors) the cartesian coordinates of
            all the atoms in the structure
        names: (list of strings) element names of all the atoms in the
            structure (used to determine whether they are hetero atoms or not)
        normal_vector: (list of 3 floats) the normal vector of the plane in its
            Hessian normal form
        point_in_plane: (list of 3 floats) the Cartesian coordinates of the
            point that lies in the plane
    """
    hydrophilic_part_indices = part_molecule(structure, normal_vector, point_in_plane)
    lipophilic_part_indices = part_molecule(
        structure, -array(normal_vector), point_in_plane
    )
    return -sum(
        [WEIGHTS_DICT[name] for name in [names[i] for i in hydrophilic_part_indices]]
    ) + sum(
        [WEIGHTS_DICT[name] for name in [names[i] for i in lipophilic_part_indices]]
    )


def pseudo_energy_complex(max_potential, point, vector, centers, potential, weights):
    """Use the complex method to determine the optimal location of the air-water interface.

    Estimate the energy the system would have if separated by a plane on whose
    left (IN the direction of the normal vector) there is an unpolar solvend
    and on whose right there is a polar solvent, like H2O. That means,
    aggregating potential unequal 0 on the right side is good whereas
    aggregating potential close to 0 is good on the left.

    Args:
        max_potential: (float) the maximum value of potential
        point: (numpy array of 3 floats) an arbitrary point on the plane that
            separates the system
        vector: (numpy array of 3 floats)  the normal vector of the plane that
            separates the system. This need not be normalized but MUST HAVE A
            NON-VANISHING NORM.
        centers: (numpy array of shape 1,3 and dtype float) Cartesian
            coordinates where the potentials are located.
        potential: (numpy array of shape 1, and dtype float) contains the
            potential values (all positive)
        weights: (numpy array of shape 1, and dtype float) if some potentials
            shall count more than others, for instance one bit should in every
            case be on a specific side, then declare something other all ones
            for this. Here, since each potential is in the center of a
            triangle, this is the area of the triangle.
    """
    left_side = dot(centers - point, vector) > 0
    right_side = invert(left_side)
    energy_left = average(potential, weights=weights * left_side)
    energy_right = -average(potential, weights=weights * right_side)
    return energy_left + energy_right


def _define_new_z_axis_gen(vector_list, new_axis, z_axis):
    """Compute the rotation matrix that can rotate the current z-axis to a new
    axis and apply this rotation matrix to a list of vectors.

    Args:
        vector_list: (list of numpy arrays of shape (3,) of dtype float) the
            vectors that shall be rotated
        new_axis: (numpy array of shape (3,) of dtype float) new z-axis
        z_axis: (numpy array of shape (3,) of dtype float) current z-axis

    Returns: 
        a generator that yields the transformed vectors (same type as before)

    See the following as an example (old state (z: old z-axis, n: new z-axis,
    vb: vector before, va: vector after))

    > z   vb
    > ↑◜◝↗
    > |φ╱
    > |╱
    > |------→ n

    The new state will be

    > z   
    > ↑
    > |
    > |
    > |------→ n
    >  ╲φ◝
    >   ╲◞
    >    ↘ va
        
    """
    rotation_axis = cross(new_axis, z_axis)
    if norm(rotation_axis) == 0:
        return (v for v in vector_list)
    [x, y, z] = rotation_axis / norm(rotation_axis)
    c = (
        1.0 * dot(new_axis, z_axis) / (norm(new_axis) * norm(z_axis))
    )  # cosine of the angle
    s = (
        1.0 * norm(cross(new_axis, z_axis)) / (norm(new_axis) * norm(z_axis))
    )  # sine of the angle
    t = 1 - c  # definition that makes the calculation much faster
    R = array(
        [
            [c + x ** 2 * t, x * y * t - z * s, x * z * t + y * s],
            [x * y * t + z * s, c + y ** 2 * t, y * z * t - x * s],
            [x * z * t - y * s, y * z * t + x * s, c + z ** 2 * t],
        ]
    )
    return (dot(R, vector) for vector in vector_list)


def get_HLB_simple(
    mol,
    thetas,
    phis,
    no_hydrogen,
    relative_angles=False,
    pseudo_energy_func=pseudo_energy_simple,
):
    """Compute the HLB value according to the simple method.

    Args:
        mol: (object of class ManipulateAggregates.aggregate.agg) the
            molecule whose HLB value shall be computed
        thetas: (list of floats) same as phis but for the theta angle of
            spherical coordinates
        phis: (list of floats) contains all the angles in the x-y-plane the
            plane should screen if relative_angles is True, this will
            correspond to the new x-y-plane after the rotation
        no_hydrogen: (bool) whether or not to irgnore bonds to hydrogens as
            cutable bonds
        relative_angles: (bool) if True, realign the z-axis to the bond vector
            prior to scanning
        pseudo_energy_func: (function like ManipulateAggregates.aggregate.hlb.pseudo_energy_simple)
            returns the pseudo energy associated with each parting plane

    Returns:
        a tuple of hlb_value, normal_vector and point_in_plane. Here, hlb_value
        (float) is the determined HLB value, normal_vector (list of 3 floats)
        is the normal vector of the plane that optimally parts the molecule
        into lipophilic and hydrophilic parts and point_in_plane (list of 3
        floats) is the point in the plane (required for Hessian normal form)

    """

    structure = array(mol.get_coordinates())
    names = mol.get_names()

    if no_hydrogen:
        # remove from the list of coordinates those that belong to hydrogen atoms since
        # the plane is not allowed to go though a bond with hydrogen
        structure_no_H = array(
            [structure[i] for i in range(0, len(names)) if names[i] != "H"]
        )
    else:
        structure_no_H = array(structure)

    # generate normal vectors allowed for the planes as vectors distributed over the
    # border of a unit sphere
    normal_vectors = [
        array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
        for theta in thetas
        for phi in phis
    ]
    normal_vectors = [
        array(v) for v in set([tuple(vector) for vector in normal_vectors])
    ]

    # the plane has to cut somewhere through a bond
    # the cutting point is taken as the middle of a bond
    # a bond is considered to exist if the two atom sare closer together than 2.55
    # Angstroms set removes duplicates but this cannot be done for lists of coordinates
    # only for tuples since lists cannot be hashed
    if relative_angles:
        # find out the coordinates of those atoms that belong to bonds and remove
        # duplicates
        bond_atoms = [
            (tuple(a), tuple(b))
            for a in structure_no_H
            for b in structure_no_H
            if norm(a - b) < 2.55
        ]
        bond_atoms = [
            [array(coord) for coord in two_coords] for two_coords in set(bond_atoms)
        ]
        # the function pseudo_energy returns an energy value per plane the lowest of
        # which corresponds to the optimum position of an interface between a
        # hydrophilic and lipophilic solvent this can be used to guess the way the
        # molecule arranges itself at such an interface
        pseudo_energy_list = []
        for bond in bond_atoms:
            bond_direction = bond[0] - bond[1]
            point = sum(bond) / 2
            for normal_vector in _define_new_z_axis_gen(
                normal_vectors, bond_direction, array([0, 0, 1])
            ):
                pseudo_energy_list.append(
                    [
                        pseudo_energy_func(structure, names, normal_vector, point),
                        normal_vector,
                        point,
                    ]
                )
    else:
        points_in_plane = [
            tuple((a + b) / 2)
            for a in structure_no_H
            for b in structure_no_H
            if norm(a - b) < 2.55
        ]
        # recast the coordinates to numpy arrays to allow easy processing
        points_in_plane = [array(coord) for coord in set(points_in_plane)]
        pseudo_energy_list = []
        for point_in_plane in points_in_plane:
            for normal_vector in _define_new_z_axis_gen(
                normal_vectors, bond_direction, array([0, 0, 1])
            ):
                pseudo_energy_list.append(
                    [
                        pseudo_energy_func(
                            structure, names, normal_vector, point_in_plane
                        ),
                        normal_vector,
                        point_in_plane,
                    ]
                )

    # find the index of the entry that contains the minimum energy
    min_index = list_min(pseudo_energy_list, key=lambda x: x[0])
    # extract the plane (as normal_vector and point_in_plane) that corresponds to the
    # minimum in energy
    normal_vector = pseudo_energy_list[min_index][1]
    point_in_plane = pseudo_energy_list[min_index][2]
    # get the indices of the atoms in the hydrophilic part
    hydrophilic_part_indices = part_molecule(structure, normal_vector, point_in_plane)
    # calculate the HLB-value after Grffin's formula as 20 times the ratio of the mass
    # in the hydrophilic part to the total mass of the molecule
    hlb_value = (
        20.0
        * sum(
            [MASSES_DICT[name] for name in [names[i] for i in hydrophilic_part_indices]]
        )
        / (sum([MASSES_DICT[name] for name in names]))
    )
    return hlb_value, normal_vector, point_in_plane


def get_HLB_complex(
    mol, thetas, phis, no_hydrogen, nr_refinements, dependence, relative_angles=False
):
    """Compute the HLB value according to the complex method.

    Args:
        mol: (object of class ManipulateAggregates.aggregate.agg) the
            molecule whose HLB value shall be computed
        thetas: (list of floats) same as phis but for the theta angle of
            spherical coordinates
        phis: (list of floats) contains all the angles in the x-y-plane the
            plane should screen if relative_angles is True, this will
            correspond to the new x-y-plane after the rotation
        no_hydrogen: (bool) whether or not to irgnore bonds to hydrogens as
            cutable bonds
        nr_refinements: (int) how many refinement steps to use during skin
            surface generation
        dependence: (string) either "dependent" (scale the potential values so
            that the highest overall [with respect to absolute value] shall be
            set to 1) or "independent" (normalize positive and negative
            potentials independently)
        phis: (list of floats) contains all the angles in the x-y-plane the
            plane should screen if relative_angles is True, this will
            correspond to the new x-y-plane after the rotation
        thetas: (list of floats) same as phis but for the theta angle of
            spherical coordinates
        relative_angles: (bool) if True, realign the z-axis to the bond vector
            prior to scanning

    Returns:
        a tuple of hlb_value, normal_vector and point_in_plane. Here, hlb_value
        (float) is the determined HLB value, normal_vector (list of 3 floats)
        is the normal vector of the plane that optimally parts the molecule
        into lipophilic and hydrophilic parts and point_in_plane (list of 3
        floats) is the point in the plane (required for Hessian normal form)
    """
    corners, face_indices, normals = mol.get_vdw_surface(
        nr_refinements=nr_refinements, shrink_factor=0.95, vdwscale=1.0
    )

    corners = array(corners, dtype=float)
    face_indices = array(face_indices, dtype=int)

    centers = average(corners[face_indices], axis=1)

    potential = array(mol.get_potential(centers))

    surfaceareas = (
        array(
            [
                norm(cross(face[1] - face[0], face[2] - face[0]))
                for face in corners[face_indices]
            ]
        )
        / 2.0
    )

    potential *= surfaceareas

    surfaceareas /= amax(fabs(surfaceareas))
    abs_potential = array(potential)
    if dependence == "dependent":
        max_potential = amax(fabs(abs_potential))
        abs_potential = fabs(abs_potential) / max_potential
    elif dependence == "independent":
        max_potential = abs(amax(abs_potential))
        min_potential = abs(amin(abs_potential))
        greater_zero = abs_potential > 0
        smaller_zero = abs_potential < 0
        abs_potential = fabs(abs_potential * greater_zero / max_potential) + fabs(
            abs_potential * smaller_zero / min_potential
        )
    else:
        raise ValueError("Dependence has to be either dependent or independent")

    coordinates = array(mol.get_coordinates())
    bondmap = mol.get_bond_map(no_hydrogen=no_hydrogen)
    if len(bondmap) == 0:
        raise RuntimeError(
            "Number of bonds through which cutting would be allowed is 0!"
        )

    # first entry is the center of the bond
    # second entry is the direction of the bond
    bonds = (
        ((coordinates[i] + coordinates[j]) * 0.5, coordinates[i] - coordinates[j])
        for i, j in bondmap
    )

    # declare theta angles
    # TODO: rework this! This is not what it's supposed to be

    # generate normal vectors allowed for the planes as vectors distributed over the border of a unit sphere
    normal_vectors = [
        (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
        for theta in thetas
        for phi in phis
    ]
    normal_vectors = [array(v) for v in set(normal_vectors)]

    # declare variables that will hold the optimum values
    min_energy = float("inf")
    normal_vector = None
    point_in_plane = None

    test_generator = (
        (point, _define_new_z_axis_gen(normal_vectors, direction, array([0, 0, 1])))
        for point, direction in bonds
    )

    max_potential = max(potential)

    for point, vectors in test_generator:
        for vector in vectors:
            energy = pseudo_energy_complex(
                max_potential, point, vector, centers, abs_potential, surfaceareas
            )
            if energy < min_energy:
                min_energy = energy
                normal_vector = vector
                point_in_plane = point

    hydromass = sum(
        mol.part_aggregate(normal_vector, point_in_plane, side="right").get_masses()
    )
    totalmass = sum(mol.get_masses())
    # calculate the HLB-value after Grffin's formula as 20 times the ratio of the mass
    # in the hydrophilic part to the total mass of the molecule
    hlb_value = 20.0 * (hydromass / totalmass)

    # part molecule
    # [i for i in range(0,len(structure)) if
    # dot(array(structure[i])-array(coordinate),array(normal_vector)) >= 0 ]

    return hlb_value, normal_vector, point_in_plane


def compute(
    mol,
    phi_step_degree=5.0,
    theta_step_degree=5.0,
    theta_range_degree=10.0,
    no_hydrogen=True,
    method="simple",
    write=False,
    nr_refinements=1,
    dependence="independent",
):
    """Wrapper function for the computation of the HLB value.

    If desired, the lipophilic and hydrophilic parts can be written to files.

    Args:
        mol: (object of class ManipulateAggregates.aggregate.agg) the
            molecule whose HLB value shall be computed
        phi_step_degree: (float) stepsize (in degrees) for phi angle for the
            plane around the bond vector
        theta_step_degree: (float) tepsize (in degrees) for theta angle for
            the plane around the bond vector
        theta_range_degree: (float) theta will vary from the negative of
            this value to the positive of this value around the bond vector
        no_hydrogen: (bool) whether or not to irgnore bonds to hydrogens as
            cutable bonds
        method: (string) "simple" or "complex" depending on which method is
            desired. See ManipulateAggregates.aggregate.hlb for more details.
        write: (bool) whether or not to write the lipophilic and hydrophilic
            parts to files (prefixed appropriately)
        nr_refinements: (int) how many refinement steps to use during skin
            surface generation
        dependence: (string) either "dependent" (scale the potential values so
            that the highest overall [with respect to absolute value] shall be
            set to 1) or "independent" (normalize positive and negative
            potentials independently). Only works if method is "complex".

    Returns:
        a tuple of hlb_value, normal_vector and point_in_plane. Here, hlb_value
        (float) is the determined HLB value, normal_vector (list of 3 floats)
        is the normal vector of the plane that optimally parts the molecule
        into lipophilic and hydrophilic parts and point_in_plane (list of 3
        floats) is the point in the plane (required for Hessian normal form)
    """

    phi_step = phi_step_degree * pi / 180.0
    theta_step = theta_step_degree * pi / 180.0
    theta_range = theta_range_degree * pi / 180.0

    thetas = concatenate(
        (
            arange(-theta_range, theta_range, theta_step),
            arange(pi + theta_range, pi - theta_range, -theta_step),
        )
    )
    phis = arange(-pi, pi, phi_step)

    # calculate the HLB-value
    # the plane that separates the hydrophilic and lipophilic parts is returned as well
    # angles are given relative to the bond
    if method == "simple":
        (hlb_value, normal_vector, point_in_plane) = get_HLB_simple(
            mol, thetas, phis, no_hydrogen, relative_angles=True
        )
    elif method == "complex":
        (hlb_value, normal_vector, point_in_plane) = get_HLB_complex(
            mol,
            thetas,
            phis,
            no_hydrogen,
            nr_refinements,
            dependence,
            relative_angles=True,
        )

    if write:
        # print both parts if desired
        filename = mol.info.get("name", "part")
        fileformat = mol.info.get("format", "xyz")
        mol.part_aggregate(normal_vector, point_in_plane, side="left").write(
            "hydrophilic_" + filename, fileformat=fileformat
        )
        mol.part_aggregate(normal_vector, point_in_plane, side="right").write(
            "lipophilic_" + filename, fileformat=fileformat
        )

    return hlb_value, normal_vector, point_in_plane
