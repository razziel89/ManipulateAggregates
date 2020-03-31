"""A useful collection of functions to read in different data files.

Supported file types are:
  - geometry: cube, molden, xyz
  - orbital data: molden
  - volumetric data: cube, dx, xyz
  - frequencies: aims, terachem
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
import re
import io
import sys
import itertools

try:
    # Python3
    import configparser as ConfigParser

    readfp = ConfigParser.ConfigParser.read_file
except ImportError:
    # Python2
    import ConfigParser

    readfp = ConfigParser.ConfigParser.readfp
from io import StringIO
import logging

logger = logging.getLogger(__name__)

try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")
try:
    from .p2p3IO import open, close, tobasestring, tounicode
except ImportError:
    logger.warning("Could not import p2p3IO")

global FLOAT_REGEX
# a regular expression matching a floating point value
FLOAT_REGEX = "^(-|)[0-9]+(\.[0-9]*){0,1}((e|E)(\+|-)[0-9]+){0,1}$"
global INT_REGEX
# a regular expression matching an value
INT_REGEX = "^(-|)[0-9]+$"


class NoOptionInConfigFileError(Exception):
    """Exception that is raised if a necessary option is not found in a config file.
    """

    pass


def read_xyz(filename):
    """Read in an xyz-file.
    
    Args:
        filename: (string) the path to the file from which data is to be read in

    Returns:
        tuple of a list of element names and a list of the Cartesian
        coordinates (lists of 3 floats)
    """
    f = open(filename, "r")

    # read the lines in the given file into the variable lines
    # and remove the trailing newline characters by using .rstrip()
    lines = numpy.array([line.rstrip() for line in f])
    # close the file descriptor
    close(f)

    # try to get the number of atoms in the molecule
    # if this does not succeed, the file is probably not a valid
    # xyz-file
    try:
        nr_atoms = int(lines[0])
    except ValueError:
        raise ValueError(
            "This is probably not a valid xyz-file since the first line does not contain an integer."
        )

    # the first two lines of an xyz file are not necessary, hence, they are removed
    # also ignore the last lines if there are more than the first line specifies
    lines = numpy.array(lines[2 : nr_atoms + 2])

    # for every line take the last three coloumns as co-ordinates for the atom
    # line.split() yields a list of the coloumns in line (whitespace separation)
    # of those, take only the last ones line.split()[1:]
    # map(float,L) gives a list whose elements are those of the list L but converted to floats
    coordinates = numpy.array([list(map(float, line.split()[1:])) for line in lines])
    # for every line take the element in the first coloumn as the name of the atom's element
    names = [line.split()[0] for line in lines]
    # if there are more entries in the file than specified in the header, ignore the additional entries
    return names, coordinates


def read_afm(filename, zscale=1):
    """Read in a file as obtained from afm measurements.

    The z-coordinate can be scaled by zscale, x- and y-coordinates will be
    multiples of 1 and all the data will be centered around 0,0

    Args:
        filename: (string) the name of the file to be read in
        zscale: (float) scale all values by this much in the z direction

    Returns:
        a list of vectors of the measured points
    """
    maximum = sys.float_info.min
    minimum = sys.float_info.max
    f = open(filename, "r")
    points = []
    (x, y) = (0.0, 0.0)
    for line in f:
        if not line.startswith("#"):
            x = 0.0
            for entry in line.rstrip().split():
                value = float(entry)
                if value > maximum:
                    maximum = value
                if value < minimum:
                    minimum = value
                p = [x, y, value]
                points.append(p)
                x = x + 1.0
            y = y + 1.0
    x -= 1
    y -= 1
    for p in points:
        p[2] = ((p[2] - minimum) / (maximum - minimum)) * zscale
        p[1] = p[1] - y / 2.0
        p[0] = p[0] - x / 2.0

    return points


def molden_positions(f, convert, elementnames=True, regex="^\[.+\]"):
    """Read the [Atoms] section of a molden file. Not for use by the user.

    Args:
        f: (file handle) a handle to the molden file
        convert: (float) every coordinate will be multiplied by this
        elementnames: (bool) whether to return the element names. Otherwise,
            the element numbers are returned.
        regex: (string) a regular expression that matches the first line that
            is no longer part of the [Atoms] section.

    Returns:
        tuple of the last line that is no longer part of the [Atoms] section
        and a list of the (possibly converted) coordinates (i.e., each a tuple
        of the element name and a list of 3 float values) in atomic units.
    """
    coords = []
    line = next(f).rstrip()
    try:
        while not re.match(regex, line):
            l = line.split()
            if elementnames:
                at_name = l[0]
            else:
                at_name = int(l[2])
            at_pos = list(map(float, l[3:6]))
            coords.append([at_name, [a * convert for a in at_pos]])
            line = next(f).rstrip()
    except (IndexError, ValueError) as e:
        raise ValueError(
            "Error on the current line: " + line + " ERROR is as follows: ", e
        )
    except StopIteration:
        pass
    return line, coords


def molden_GTO(
    f,
    GTO_coefficients=False,
    nr_primitives=True,
    regex="^\[.+\]",
    int_regex=INT_REGEX,
    float_regex=FLOAT_REGEX,
    orbital_regex="^(s|p|d|f|g|S|P|D|F|G)$",
):
    """Read GTO section of a molden file. Not for use by the user.

    Args:
        f: (file handle) a handle to the molden file
        GTO_coefficients: see ManipulateAggregates.collection.read.read_molden
        nr_primitives: same as GTO_nr_primitives for ManipulateAggregates.collection.read.read_molden
        regex: (string) a regular expression declaring the start of any
            section. The section is considered to end either and EOF or when
            this regex matches the current line.
        int_regex: (string) a regex matching an integer number
        float_regex: (string) a regex matching a float number
        orbital_regex: (string) a regex matching all valid orbital names

    Returns:
        tuple of the last line that is no longer part of the [GTO] section and
        a list of the Gaussian type orbitals (GTO). A GTO is a list of at least
        one integer (the meaning of this depends on nr_primitives). If
        GTO_coefficients is True, also a list of shells will be returned. A
        shell is a list of one letter (the shell's type), one float (a number
        required for backwards compatibility) and one int (the number of
        primitives in the shell) and a list of primitives. A primitive is a
        list of 2 floats.  First the exponential factor and second the
        prefactor.
    """
    dict_orb_numbers = {"s": 1, "p": 3, "d": 6, "f": 10, "g": 15}
    gto = []
    at_gto = 0
    if GTO_coefficients:
        at_gto_coeff = []
        shell = []
    at_nr = 0
    count = 0
    elements = 0
    line = next(f).rstrip()
    try:
        while not re.match(regex, line):
            l = line.split()
            if len(l) == 2 and re.match(int_regex, l[0]) and re.match(int_regex, l[1]):
                if GTO_coefficients or nr_primitives:
                    if not count == elements:
                        raise ValueError(
                            "Atom "
                            + str(at_nr)
                            + " needs "
                            + str(elements)
                            + " elements for basis function but only "
                            + str(count)
                            + " are available. "
                            + line
                        )
                if at_gto > 0 and at_nr > 0:
                    if GTO_coefficients:
                        at_gto_coeff.append([orbitaltype, stretch, elements, shell])
                        shell = []
                        gto.append([at_gto, at_gto_coeff])
                        at_gto_coeff = []
                    else:
                        gto.append(at_gto)
                at_gto = 0
                at_nr += 1
            elif (
                len(l) == 3
                and re.match(orbital_regex, l[0])
                and re.match(int_regex, l[1])
                and re.match(float_regex, l[2])
            ):
                if GTO_coefficients:
                    if len(shell) > 0:
                        at_gto_coeff.append([orbitaltype, stretch, elements, shell])
                    shell = []
                orbitaltype = l[0]
                stretch = float(l[2])
                if GTO_coefficients or nr_primitives:
                    elements = int(l[1])
                    count = 0
                if nr_primitives:
                    at_gto += dict_orb_numbers[l[0].lower()] * elements
                else:
                    at_gto += dict_orb_numbers[l[0].lower()]
            elif (
                len(l) == 2
                and re.match(float_regex, l[0])
                and re.match(float_regex, l[1])
            ):
                if nr_primitives or GTO_coefficients:
                    count += 1
                if GTO_coefficients:
                    tempshell = list(map(float, [l[0], l[1]]))
                    shell.append(tempshell)
            elif re.match("^\s*$", line):
                pass
            line = next(f).rstrip()
    except (IndexError, ValueError) as e:
        raise ValueError(
            "Error on the current line: " + line + " ERROR is as follows: ", e
        )
    except StopIteration:
        pass
    if at_nr > 0:
        if GTO_coefficients:
            at_gto_coeff.append([orbitaltype, stretch, elements, shell])
            sheel = []
            gto.append([at_gto, at_gto_coeff])
        else:
            gto.append(at_gto)
    return line, gto


def molden_MO(
    f,
    MO_coefficients=False,
    regex="^\[.+\]",
    int_regex=INT_REGEX,
    float_regex=FLOAT_REGEX,
):
    """Read MO section of a molden file. Not for use by the user.

    Args:
        f: (file handle) a handle to the molden file
        MO_coefficients: see ManipulateAggregates.collection.read.read_molden
        regex: (string) a regular expression declaring the start of any
            section. The section is considered to end either and EOF or when
            this regex matches the current line.
        int_regex: (string) a regex matching an integer number
        float_regex: (string) a regex matching a float number

    Returns:
        tuple of the last line that is no longer part of the [MO] section and
        a list of the molecular orbitals (MO). A MO is a list of at least
        one float (the MO's energy), one string ('alpha' or 'beta' depending on
        the spin) and one float (the MO's occupation number). If
        MO_coefficients, also a list of floats will be returned (the linear
        combination factors that, together with the shells, make up the MO).
    """
    mo = []
    orbital = []
    energy = None
    # is_alpha==True means it is alpha, False means it's beta
    is_alpha = None
    occupation = None
    line = next(f).rstrip()
    try:
        while not re.match(regex, line):
            l = line.split()
            if len(l) == 2 and re.match("^.+=$", l[0]):
                if MO_coefficients:
                    if len(orbital) > 0:
                        if not (
                            occupation is None or energy is None or is_alpha is None
                        ):
                            mo.append(
                                [
                                    energy,
                                    "alpha" if is_alpha else "beta",
                                    occupation,
                                    orbital,
                                ]
                            )
                            energy = None
                            is_alpha = None
                            occupation = None
                        else:
                            raise ValueError(
                                "At least one orbital does not specify all of Ene=, Spin= and Occup="
                            )
                elif (
                    not occupation == None
                    and not energy == None
                    and not is_alpha == None
                ):
                    mo.append([energy, "alpha" if is_alpha else "beta", occupation])
                    energy = None
                    is_alpha = None
                    occupation = None
                if MO_coefficients:
                    orbital = []
                    count = 0
                if l[0] == "Ene=":
                    energy = float(l[1])
                elif l[0] == "Spin=":
                    if re.match("^[Aa][Ll][Pp][Hh][Aa]$", l[1]):
                        is_alpha = True
                    elif re.match("^[Bb][Ee][Tt][Aa]$", l[1]):
                        is_alpha = False
                    else:
                        raise ValueError(
                            "Spin name has to be either alpha or beta but it is " + l[1]
                        )
                elif l[0] == "Occup=":
                    occupation = float(l[1])
            elif (
                len(l) == 2
                and re.match(int_regex, l[0])
                and re.match(float_regex, l[1])
            ):
                if MO_coefficients:
                    newcount = int(l[0])
                    # if contributions have been left out from the molden file, fill them with zeroes
                    for i in range(count + 1, newcount):
                        orbital.append(0.0)
                    count = newcount
                    orbital.append(float(l[1]))
            elif re.match("^\s*$", line):
                pass
            line = next(f).rstrip()
    except (IndexError, ValueError, ValueError) as e:
        raise ValueError(
            "Error on the current line: " + line + " ERROR is as follows: ", e
        )
    except StopIteration:
        pass
    if MO_coefficients:
        if len(orbital) > 0:
            if not (occupation is None or energy is None or is_alpha is None):
                mo.append(
                    [energy, "alpha" if is_alpha else "beta", occupation, orbital]
                )
            else:
                raise ValueError(
                    "At least one orbital does not specify all of Ene=, Spin= and Occup="
                )
    elif not (occupation is None or energy is None or is_alpha is None):
        mo.append([energy, "alpha" if is_alpha else "beta", occupation])

    if MO_coefficients:
        max_nr_primitives = max([len(o[-1]) for o in mo])
        # set all orbitals that have been filled to the same maximum range
        for i in range(len(mo)):
            mo[i][-1] += [0.0] * (max_nr_primitives - len(mo[i][-1]))
    return line, mo


def read_molden(
    filename,
    positions=True,
    elementnames=True,
    GTO=True,
    GTO_coefficients=False,
    GTO_nr_primitives=False,
    MO=True,
    MO_coefficients=False,
):
    """Read in a Molden file accordig to some flags that are set.

    Beware: only Cartesian coordinates are supported so far. Will return a
    dictionary with approrpiately named entries. Only Cartersian Gaussian type
    orbitals are supported as of now.

    All units will be converted to Bohr when reading in the file. This
    function assumes that all exponents in the GTO section are given in Bohr^-1
    (which they usually are).

    Args:
        filename: (string) the path to the file from which data is to be read in
        positions: (bool) whether or not atomic coordinates shall be read in.
        elementnames: (bool) whether or not elements shall be identified by
            their names (True) or by their element numbers (False)
        GTO: (bool) whether or not to read in the GTO section
        GTO_coefficients: (bool) whether or not to read in all
            GTO-coefficients. If False, only the number of primitives or
            shells are counted and returned (influenced by GTO_nr_primitives)
        GTO_nr_primitives: (bool) whether you want to count the number of
            primitives (if True) or shells (if False)
        MO: (bool) whether or not to read in the MO-section. You will get 2 lists
            one for each spin containing energy and occupation
        MO_coefficients: (bool) whether or not all MO-coefficients shall be
            read in If False, only auxilliary information line energies and
            occupations and spins are read in.

    Returns:
        a dictionary with apprropriately named entries. Keys are "GTO"
        (via ManipulateAggregates.collection.read.molden_GTO), "positions"
        (via ManipulateAggregates.collection.read.molden_positions) and "MO"
        (via ManipulateAggregates.collection.read.molden_MO). Please see the
        documentation of the return values of these functions for the format of
        the returned arguments. Note that only the second argument of the
        tuples these functions returned is part of the dictionary that this
        function returns.
    """
    # regular expression that matches the beginning of any secion
    regex = "^\[.+\]"
    # regular expression that matches an integer
    global INT_REGEX
    int_regex = INT_REGEX
    # regular expression that matches a float
    global FLOAT_REGEX
    float_regex = FLOAT_REGEX
    # regular expression that matches any of the allowed orbital names
    orbital_regex = "^(s|p|d|f|g|S|P|D|F|G)$"
    # open the file
    f = open(filename, "r")
    # initialize dictionary that will hold the results
    result = {}
    # check format of first line
    if not re.match("^\[Molden Format\]\s*$", next(f).rstrip()):
        raise ValueError("The first line in a Molden file has to be '[Molden Format]'")
    nr_sections = [positions, GTO, MO].count(True)
    sec = nr_sections - 1
    important_sections_regex = ""
    if positions:
        important_sections_regex += "^\[Atoms\]\s+([Aa][Nn][Gg][Ss]|[Aa][Uu])\s*$"
        if sec > 0:
            important_sections_regex += "|"
            sec -= 1
    if GTO:
        important_sections_regex += "^\[GTO\]"
        if sec > 0:
            important_sections_regex += "|"
            sec -= 1
    if MO:
        important_sections_regex += "^\[MO\]"
        if sec > 0:
            important_sections_regex += "|"
            sec -= 1
    sec = 0
    try:
        # skip to the next important section
        line = next(f).rstrip()
        while not re.match(important_sections_regex, line):
            line = next(f).rstrip()
        # read in only the requested number of sections
        while sec < nr_sections:
            # read in requested sections
            if (
                re.match("^\[Atoms\]\s+([Aa][Nn][Gg][Ss]|[Aa][Uu])\s*$", line)
                and positions
            ):
                # Atoms section
                sec += 1
                # determine whether the atomic coordinates are given in angstroms or bohrs
                if re.match("^\[Atoms\]\s+([Aa][Nn][Gg][Ss])\s*$", line):
                    convert = 1.0 / 0.52918
                else:
                    convert = 1.0
                # read in the section
                line, result["positions"] = molden_positions(
                    f, convert, elementnames=elementnames, regex=regex
                )
            elif re.match("^\[GTO\]", line) and GTO:
                # GTO section
                sec += 1
                # read in the section
                line, result["GTO"] = molden_GTO(
                    f,
                    GTO_coefficients=GTO_coefficients,
                    regex=regex,
                    nr_primitives=GTO_nr_primitives,
                    int_regex=int_regex,
                    float_regex=float_regex,
                    orbital_regex=orbital_regex,
                )
            elif re.match("^\[MO\]", line) and MO:
                # MO section
                # GTO section
                sec += 1
                # read in the section
                line, result["MO"] = molden_MO(
                    f,
                    MO_coefficients=MO_coefficients,
                    regex=regex,
                    int_regex=int_regex,
                    float_regex=float_regex,
                )

            # break loop prematurely if enough sections have been read in
            if re.match(regex, line):
                if sec >= nr_sections:
                    break
            # skip to next significant section if it has not yet been reached
            if not re.match(important_sections_regex, line):
                while not re.match(important_sections_regex, line):
                    line = next(f).rstrip()
    except StopIteration:
        if sec < nr_sections:
            # if end of file reached before the requested sections could be read in
            raise ValueError(
                "You requested "
                + str(nr_sections)
                + " sections but only "
                + str(sec)
                + " were found as end of file was reached."
            )
    # adjust the number of entries per MO so that there is a coefficient for each primitive
    if MO and GTO:
        if GTO_nr_primitives:
            nr_primitives = sum(gto[0] for gto in result["GTO"])
        else:
            NR_PRIMS = {
                "s": 1,
                "p": 3,
                "d": 6,
                "f": 10,
            }
            if GTO_coefficients:
                nr_primitives = sum(
                    NR_PRIMS[shell[0]] for gto in result["GTO"] for shell in gto[1]
                )
            else:
                # since the actual number of primitives cannot be known unless GTO_coefficients is True,
                # the actual length cannot be properly adjusted (but might not be necessary)
                nr_primitives = 0
        for i in range(len(result["MO"])):
            len_diff = nr_primitives - len(result["MO"][i][-1])
            # if len_diff<0, nothing will be added
            result["MO"][i][-1] += [0.0] * len_diff
    # do some sanity checks on the results
    if positions and len(result["positions"]) == 0:
        raise ValueError(
            "Atomic coordinates requested but whole file read in without finding the secion."
        )
    if GTO and len(result["GTO"]) == 0:
        raise ValueError(
            "GTO section requested but whole file read in without finding the secion."
        )
    if MO and len(result["MO"]) == 0:
        raise ValueError(
            "MO section requested but whole file read in without finding the secion."
        )
    close(f)
    return result


def _renormalize_whole(vec, norm=1.0):
    return vec * (norm / numpy.linalg.norm(vec))


def _renormalize_individual(vec, norm=1.0):
    return vec * (norm / numpy.max(numpy.linalg.norm(vec, axis=1)))


def _renormalize_none(vec, norm=1.0):
    return vec * norm


def read_aims_frequencies(fname, mode=1, amplitude_factor=1, normalize="individual"):
    """Read in a frequencies file in AIMS format.

    Displacement vectors can be normalized.

    Args:
        fname: (string) filename
        mode: (int) mode to be read in (starting at 1)
        amplitude_factor: (float) the new norm of the vector according to the
            value of normalize
        normalize: possible values: whole, individual, none. whole
            (normalize the entire displacement vector), individual (for the
            atom with the largest total displacement vector, normalize this to
            amplitude_factor and adjust all the others accordingly.), none
            (do not perform normalization, only scale by amplitude_factor)

    Returns:
        a tuple of 2 lists. The first list (of floats) contains the frequencies
        associated with the modes. The second list (actually a numpy array with
        a shape of N,3 [N being the number of atoms]) contains the
        displacements for each atom for the corresponding mode.
    """
    # read the lines in the given file into the variable lines
    # and remove the trailing newline characters by using .rstrip()
    f = open(fname, "r")
    line1 = next(f).rstrip()

    nr_atoms = int(line1.split()[0])

    nr_deg_of_freedom = 3 * nr_atoms
    nr_normalmodes = int(line1.split()[1])

    displacement = numpy.zeros(nr_deg_of_freedom)

    mode_count = 0
    while mode_count < mode:
        line = next(f).rstrip()
        freqs = map(float, line.split())
        for value in freqs:
            mode_count += 1
            if mode_count == mode:
                frequency = value
    while mode_count < nr_normalmodes:
        mode_count += len(next(f).rstrip().split())

    coord = 0
    while coord < mode * nr_deg_of_freedom:
        line = next(f).rstrip()
        disp = map(float, line.split())
        for value in disp:
            if (
                coord >= (mode - 1) * nr_deg_of_freedom
                and coord < mode * nr_deg_of_freedom
            ):
                displacement[coord % nr_deg_of_freedom] = value
            coord += 1
    close(f)
    displacement.shape = (nr_atoms, 3)
    if normalize == "individual":
        displacement = _renormalize_individual(displacement, amplitude_factor)
    elif normalize == "whole":
        displacement = _renormalize_whole(displacement, amplitude_factor)
    else:
        displacement = _renormalize_none(displacement, amplitude_factor)

    return frequency, displacement


def read_terachem_frequencies(
    fname, mode=1, amplitude_factor=1, normalize="individual"
):
    """Read in a frequencies file in TeraChem format.

    Args:
        fname: (string) filename
        mode: (int) mode to be read in (starting at 1). If 0 given, all modes
            are read in (normalization not supported for mode == 0)
        amplitude_factor: (float) the new norm of the vector according to the value
            of normalize
        normalize: possible values: whole, individual, none. whole
            (normalize the entire displacement vector), individual (for the
            atom with the largest total displacement vector, normalize this to
            amplitude_factor and adjust all the others accordingly.), none
            (do not perform normalization, only scale by amplitude_factor).
            NOTE that normalization does not work for mode == 0

    Returns:
        a tuple of 2 lists. The first list (of floats) contains the frequencies
        associated with the modes. The second list (actually a numpy array with
        a shape of N,3 [N being the number of atoms]) contains the
        displacements for each atom for the corresponding mode. If mode is
        0, each element of the tuple is a list of the aforementioned.
    """
    f = open(fname, "r")
    # read the lines in the given file into the variable lines
    # and remove the trailing newline characters by using .rstrip()
    line1 = next(f).rstrip()

    nr_atoms = int(line1.split()[2])

    nr_deg_of_freedom = 3 * nr_atoms

    # number of normalmodes is in second line
    line1 = next(f).rstrip()
    nr_normalmodes = int(line1.split()[3])

    # the third line does not contain any useful information
    next(f)

    if mode == 0:
        displacement = numpy.zeros((nr_deg_of_freedom, nr_normalmodes), dtype=float)
        freq = numpy.zeros(nr_normalmodes, dtype=float)

        maxmode = 0
        while maxmode < nr_normalmodes:
            # skip this line
            next(f)
            entries = next(f).rstrip().split()
            freq[maxmode : maxmode + len(entries)] = list(map(float, entries))
            # skip this line
            next(f)
            for i in range(nr_atoms):
                # this line contains the number of the atom and the x-displacement for
                # the first len(entries)-1 modes
                entries = next(f).rstrip().split()
                at = int(entries[0])
                displacement[(at - 1) * 3 + 0][
                    maxmode : maxmode + len(entries) - 1
                ] = list(map(float, entries[1:]))
                # this line contains the y-displacement for the first len(entries) modes
                entries = next(f).rstrip().split()
                displacement[(at - 1) * 3 + 1][maxmode : maxmode + len(entries)] = list(
                    map(float, entries)
                )
                # this line contains the z-displacement for the first len(entries) modes
                entries = next(f).rstrip().split()
                displacement[(at - 1) * 3 + 2][maxmode : maxmode + len(entries)] = list(
                    map(float, entries)
                )
            maxmode += len(entries)
            next(f)
        norm = numpy.linalg.norm(displacement, axis=0)
        norm.shape = (1, nr_normalmodes)
        displacement /= norm
    else:
        displacement = numpy.zeros(nr_deg_of_freedom, dtype=float)
        disp_count = 0
        # mode counting starts at 0
        mode -= 1
        maxmode = -1
        while maxmode < mode:
            modes = list(map(int, next(f).rstrip().split()))
            freqs = list(map(float, next(f).rstrip().split()))
            # the next line does not contain any useful information
            next(f)
            maxmode = max(modes)
            if maxmode < mode:
                for i in range(nr_atoms * 3 + 1):
                    next(f)
            else:
                index = modes.index(mode)
                freq = freqs[index]
                for at in range(nr_atoms):
                    # treat first line per atom
                    entries = next(f).rstrip().split()
                    displacement[disp_count] = list(map(float, entries[1:])[index])
                    disp_count += 1
                    # treat second line per atom
                    entries = next(f).rstrip().split()
                    displacement[disp_count] = list(map(float, entries)[index])
                    disp_count += 1
                    # treat third line per atom
                    entries = next(f).rstrip().split()
                    displacement[disp_count] = list(map(float, entries)[index])
                    disp_count += 1
        displacement.shape = (nr_atoms, 3)
        if normalize == "individual":
            displacement = _renormalize_individual(displacement, amplitude_factor)
        elif normalize == "whole":
            displacement = _renormalize_whole(displacement, amplitude_factor)
        else:
            displacement = _renormalize_none(displacement, amplitude_factor)

    close(f)

    return freq, displacement


def read_charges_simple(filename, compare_elements=False, molecule=None):
    """Read in an xyz-file where each line of Cartesian coordinates is followed by a charge.
    
    Args:
        filename: the path to the file from which data is to be read in
        compare_elements: (bool) if this is True, molecule must be of type
            ManipulateAggregates.aggregate.agg. Then a sanity check will be
            performed where the element names from the molecule object are compared
            to those from the given file. Furthermore, the coordinates are taken
            not from the file but from the molecule object.
        molecule: (ManipulateAggregates.aggregate.agg) the molecule object
            to compare against if compare_elements is True.

    Returns:
        the first argument to be returned will be the position of the partial
        charges (as lists of 3 floats) and the second will be a list of the
        charges (as floats).
    """
    f = open(filename, "r")

    # read the lines in the given file into the variable lines
    # and remove the trailing newline characters by using .rstrip()
    lines = numpy.array([line.rstrip().split() for line in f])
    # close the file descriptor
    close(f)

    # try to get the number of atoms in the molecule
    # if this does not succeed, the file is probably not a valid
    # xyz-file
    try:
        nr_atoms = int(lines[0][0])
    except ValueError:
        raise ValueError(
            "This is probably not a valid xyz-file since the first line does not contain an integer."
        )

    # the first two lines of an xyz file are not necessary, hence, they are removed
    # also ignore the last lines if there are more than the first line specifies
    lines = numpy.array(lines[2 : nr_atoms + 2])

    if compare_elements:
        elements_molecule = [e for e in sorted(molecule.get_names())]
        elements_file = [line[0] for line in sorted(lines)]
        if elements_molecule == elements_file:
            charges = numpy.array([float(line[4]) for line in lines])
            coordinates = numpy.array(molecule.get_coordinates())
        else:
            raise ValueError(
                "Molecule read from the charge file and the given molecule object do not contain the same elements."
            )
    else:
        try:
            coordinates = numpy.array([list(map(float, line[1:4])) for line in lines])
            charges = numpy.array([float(line[4]) for line in lines])
        except IndexError:
            raise IndexError(
                "Not enough coloumns! There need to be 4 in every line but the first to. Element name, x,y,z coordinates, charge."
            )
        except ValueError:
            raise ValueError(
                "At least one value on one of the lines is no valid float."
            )
    return coordinates, charges


def _list_equiv(l1, l2):
    """Check if two lists are equivalent by comparing them element-wise. If
    one list is exhausted and they were equivalent so far, consider them
    equivalent.
    """
    for e1, e2 in zip(l1, l2):
        if not e1 == e2:
            return False
    return True


def read_charges_dx(
    filename,
    add_nuclear_charges=False,
    molecule=None,
    unit_conversion=1.0,
    invert_charge_data=True,
    rescale_charges=True,
    total_charge=0,
    nr_return=None,
    density=False,
    header_dict=None,
):
    """Read in a DX file.
    
    Args:
        filename: (string) the name of the dx file
        add_nuclear_charges: (bool) if True, nuclear charges and coordinates will
            be the first entries in the returned lists, needs molecule to be declared
        molecule: (object of class ManipulateAggregates.aggregate.agg)
            the get_coordinates() and get_charges() methods are used to get the
            coordinates and values of nuclear point charges.
        unit_conversion: (float) a value by which to scale all coordinates
        invert_charge_data: (bool) if False, do not invert the volumetric
            charge data. Nuclear charges are always positive and volumetric
            data is taken inverted by default.
        rescale_charges: (bool) the sum of atomic charges and the sum of
            nuclear charges have to match. If this is True and the charges
            don't match, rescale all volumetric data linearly so that they do.
            Only makes sense if add_nuclear_charges == True
        total_charge: (float) the total charge of the molecule to properly
            rescale the electronic charges (if requested)
        nr_return: (list or None) if a variable of type list is given, append
            to it the number of atoms and the number of volumetric entries
        density: (bool) if True, return the density at the center of the voxel
            instead of the product of the density and the voxel's volume
        header_dict: (dictionary or None) if given a dictionary, add to it the
            entries "counts_xyz", "org_xyz", "delta_x", "delta_y" and "delta_z"
            to allow easy re-printing of the dx-file

    Returns:
        tuple of 2 lists. The first list contains 3-element lists of floats
        that are the Cartesian coordinates of the volumetric data values. The
        second list is a flat list of floats that contains the associated
        volumetric data.
    """
    dxdata = read_dx(
        filename,
        unit_conversion=unit_conversion,
        invert_charge_data=invert_charge_data,
        density=density,
        header_dict=header_dict,
        grid=True,
        data=True,
    )
    charges = numpy.array(dxdata["data"])
    coordinates = numpy.array(dxdata["grid"])

    if add_nuclear_charges:
        if not molecule is None:
            try:
                # since no data about the nuclei is saved in the dx file,
                # get it from the molecule object
                nuc_charges = numpy.array(molecule.get_charges())
                nuc_coordinates = numpy.array(molecule.get_coordinates())

            except AttributeError as e:
                raise AttributeError(
                    "Given molecule object does not define methods get_charges or get_coordinates.",
                    e,
                )
        else:
            raise ValueError("Cannot add nuclear charges without molecule argument.")

        if rescale_charges:
            # This equatiom makes it so that the sum over the electronic charges (which is negative)
            # plus the sum of the nuclear charges (which is positive) equals the total charge.
            charges *= (total_charge - numpy.sum(nuc_charges)) / numpy.sum(charges)

        numpy.insert(coordinates, 0, nuc_coordinates, axis=0)
        numpy.insert(charges, 0, nuc_charges, axis=0)

    return coordinates, charges


def gziplines(fname):
    """Return a generator used to read a gzipped file line-wise.

    Args:
        fname: (string) the filename

    Returns:
        a generator
    """
    from subprocess import Popen, PIPE

    f = Popen(["zcat", fname], stdout=PIPE, bufsize=4096)
    for line in f.stdout:
        yield tobasestring(line)
    if f.wait() != 0:
        raise ValueError("Error when reading in gzipped file%s" % (fname))


_filetype_regexes = {
    "gzipped": re.compile(r"^gzip compressed data\b", re.IGNORECASE),
    "text": re.compile(r"^(utf-[0-9]+|ascii).* text\b", re.IGNORECASE),
}


def is_gzipped(fname):
    """Determine whether or not a file is gzipped.

    This uses the "file" programme. Sometimes, this programme does not
    correctly determine that a file is a gzipped plain text file. In such a
    case, this function tries to read from the file assuming it is gzipped. If
    that succeeds, True is returned.

    Args:
        fname: (string) the filename

    Returns:
        whether or not the file is gzipped.
    """
    from subprocess import Popen, PIPE

    f = Popen(["file", "-n", "-b", fname], stdout=PIPE, bufsize=4096)
    gzipped = False
    lines = 0
    for dline in f.stdout:
        line = tobasestring(dline)
        lines += 1
        if re.match(_filetype_regexes["gzipped"], line) is not None:
            gzipped = True
        elif re.match(_filetype_regexes["text"], line) is not None:
            gzipped = False
        else:
            print(
                "WARNING: Could not determine whether file is gzipped or not from magic bytes, will try differently",
                file=sys.stderr,
            )
            print("         File is: %s" % (fname), file=sys.stderr)
            import gzip

            handle = gzip.open(fname, "rb")
            try:
                next(handle)
                print("         File determined to be gzipped.", file=sys.stderr)
                gzipped = True
            except (IOError, OSError) as e:
                print(
                    "         File probably not gzipped, expect reading it to fail though, IOError was:",
                    file=sys.stderr,
                )
                print(e, file=sys.stderr)
                gzipped = False
            finally:
                handle.close()
    if lines != 1:
        raise IOError("The 'file' command gave more than one line of output.")
    f.wait()
    return gzipped


def read_dx(
    filename,
    unit_conversion=1.0,
    invert_charge_data=False,
    density=True,
    header_dict=None,
    grid=True,
    data=True,
    silent=False,
    gzipped=False,
    comments=False,
):
    """Read in a DX file.
    
    Args:
        filename: (string) the file name
        unit_conversion: (float) give a value by which to scale all coordinates
        invert_charge_data: (bool) if False, do not invert the volumetric
            charge data
        density: (bool) if True, return the density at the center of the voxel
            instead of the product of the density and the voxel's volume
        header_dict: (empty dictionary or None) if given a dictionary, add to
            it the entries "counts_xyz", "org_xyz", "delta_x", "delta_y" and
            "delta_z" to allow easy re-printing of the dx-file.
        grid: (bool) whether or not to return the grid.
        data: (bool) whether or not to return the data.
        silent: (bool) whether or not to utter warnings such as about a missing footer
        gzipped: (bool) treat the file as a gzipped one. If it cannot be
            treated as such, default to non-gzipped mode.
        comments: (bool) whether or not to return comments saved in the dx file

    Returns:
        a dictionary with the entries "grid" (if grid == True) and "data"
        (if data == True) and "comments" at start of file (if comments == True)
    """
    result = {}
    if gzipped and is_gzipped(filename):
        f = gziplines(filename)
    else:
        f = open(filename, "r")
    if header_dict is not None:
        import copy
    l = next(f)
    if comments:
        # save comment lines at the beginning
        result["comments"] = []
        while l.startswith("#"):
            result["comments"].append(l.rstrip())
            l = next(f)
    else:
        # ignore comment lines at the beginning
        while l.startswith("#"):
            l = next(f)
    # according to the VMD mailing list, the ordering is: z fast, y medium, and x slow
    # this translates to: z inner, y middle, and x outer
    # Source: http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/21526.html
    match_dict = {"outer": 0, "inner": 2, "middle": 1}
    xyz_dict = {"x": 0, "y": 1, "z": 2}
    # preallocate stuff
    axes = [None] * 3
    nrs = [None] * 3
    # Comment on units: the grid (if requested) will have the units used in the file itself multiplied
    #                  by unit_conversion
    # read in the header
    # line stating the number of gridpositions
    line = l.rstrip().split()
    if _list_equiv(line, ["object", "1", "class", "gridpositions", "counts"]):
        try:
            nrs = list(map(int, line[5:]))
        except ValueError as e:
            raise ValueError(
                "First non-comment line in DX file does not end on three integers. Line: %s"
                % (line),
                e,
            )
    else:
        raise ValueError(
            "First non-comment line in DX file must be 'object 1 class gridpositions counts nx ny nz' Line: %s"
            % (line)
        )
    if header_dict is not None:
        header_dict["counts_xyz"] = copy.copy(nrs)
    # line with origin
    line = next(f).rstrip().split()
    if line[0] == "origin":
        try:
            origin = numpy.array(list(map(float, line[1:]))) * unit_conversion
        except ValueError as e:
            raise ValueError(
                "Second non-comment line in DX file does not end on three floats. Line: %s"
                % (line),
                e,
            )
    else:
        raise ValueError(
            "Second non-comment line in DX file must be 'origin ox oy oz' Line: %s"
            % (line)
        )
    if header_dict is not None:
        header_dict["org_xyz"] = copy.copy(list(origin))

    # the three lines with the voxel sizes
    for c in range(3):
        line = next(f).rstrip().split()
        if line[0] == "delta":
            try:
                axes[c] = list(map(float, line[1:4]))
            except ValueError as e:
                raise ValueError(
                    "One of third to fifth non-comment lines in DX file does not end on three floats. Line: %s"
                    % (line),
                    e,
                )
        else:
            raise ValueError(
                "Third to fifth non-comment lines must be 'delta dx dy dz' Line: %s"
                % (line)
            )
    if header_dict is not None:
        header_dict["delta_x"] = copy.copy(axes[0])
        header_dict["delta_y"] = copy.copy(axes[1])
        header_dict["delta_z"] = copy.copy(axes[2])
    # next line, which somewhat of a duplicate, the integers are ignored here
    line = next(f).rstrip().split()
    if _list_equiv(line, ["object", "2", "class", "gridconnections", "counts"]):
        try:
            list(map(int, line[5:]))
        except ValueError as e:
            raise ValueError(
                "Sixth non-comment line in DX file does not end on three integers. Line: %s"
                % (line),
                e,
            )
    else:
        raise ValueError(
            "Sixth non-comment line in DX file must be 'object 2 class gridconnections counts nx ny nz' Line: %s"
            % (line)
        )
    # next line contains some test data to check whether format is correct
    line = next(f).rstrip().split()
    if _list_equiv(
        line, ["object", "3", "class", "array", "type", "double", "rank", "0", "items"]
    ) and _list_equiv(line[10:], ["data", "follows"]):
        try:
            if not nrs[0] * nrs[1] * nrs[2] == int(line[9]):
                raise ValueError(
                    "Seventh non-comment line in DX file does not contain the correct number of volumetric data elements. Line: %s"
                    % (line)
                )
        except ValueError as e:
            raise ValueError(
                "Seventh non-comment line in DX file does not contain the total number of volumetric data elements Line: %s"
                % (line),
                e,
            )
    else:
        raise ValueError(
            "Seventh non-comment line in DX file must be 'object 3 class array type double rank 0 items nx*ny*nz data follows' Line: %s"
            % (line)
        )
    # convert to numpy arrays
    axes = numpy.array(axes)
    nrs = numpy.array(nrs)
    # this is the volume of one voxel
    # which will be used to convert charge density to charge
    # so it will seem as if there were a point charge at the center
    # of the voxel containing the whole charge inside that voxel
    if density:
        volume = 1.0
    else:
        volume = numpy.linalg.det(unit_conversion * axes.transpose())
    # read in volumetric data if requested
    if data:
        charges = numpy.zeros((nrs[0] * nrs[1] * nrs[2]), dtype=float)
        count = 0
        for l in f:
            if l.startswith(
                ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "+", "-")
            ):
                for e in map(float, l.rstrip().split()):
                    charges[count] = e * volume
                    count += 1
            else:
                break
        if invert_charge_data:
            charges *= -1
        result["data"] = charges
    else:
        for l in f:
            if not l.startswith(
                ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "+", "-")
            ):
                break
    try:
        # read in footer which is assumed to be of a certain format
        # this might break if a programme changes data assignments
        # in l is the first line of the footer
        while l.startswith(
            ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "+", "-")
        ):
            l = next(f)
        l = l.rstrip().split()
        if not _list_equiv(l, ["attribute", '"dep"', "string", '"positions"']):
            raise ValueError(
                'First line of footer must be attribute \'"dep" string "positions"\' Line: %s'
                % (l)
            )
        l = next(f).rstrip().split()
        if not _list_equiv(
            l,
            [
                "object",
                '"regular',
                "positions",
                "regular",
                'connections"',
                "class",
                "field",
            ],
        ):
            raise ValueError(
                "Second line of footer must be attribute 'object \"regular positions regular connections\" class field' Line: %s"
                % (l)
            )
        l = next(f).rstrip().split()
        if not _list_equiv(l, ["component", '"positions"', "value", "1"]):
            raise ValueError(
                "Third line of footer must be attribute 'component \"positions\" value 1' Line: %s"
                % (l)
            )
        l = next(f).rstrip().split()
        if not _list_equiv(l, ["component", '"connections"', "value", "2"]):
            raise ValueError(
                "Third line of footer must be attribute 'component \"connections\" value 2' Line: %s"
                % (l)
            )
        l = next(f).rstrip().split()
        if not _list_equiv(l, ["component", '"data"', "value", "3"]):
            raise ValueError(
                "Third line of footer must be attribute 'component \"data\" value 3' Line: %s"
                % (l)
            )
    except StopIteration:
        if not silent:
            print("WARNING: no footer present in dx-file", file=sys.stderr)
    # preallocate
    if grid:
        coordinates = numpy.zeros((nrs[0] * nrs[1] * nrs[2], 3), dtype=float)
        axes_rearranged = numpy.zeros(axes.shape, dtype=float)
        axes_rearranged[0] = axes[match_dict["outer"]]
        axes_rearranged[1] = axes[match_dict["middle"]]
        axes_rearranged[2] = axes[match_dict["inner"]]
        nrs_rearranged = numpy.zeros(nrs.shape, dtype=int)
        nrs_rearranged[0] = nrs[match_dict["outer"]]
        nrs_rearranged[1] = nrs[match_dict["middle"]]
        nrs_rearranged[2] = nrs[match_dict["inner"]]
        # this is a fast variant to do the rearranging of the data. see
        # read_charges_cube for the other variants
        result["grid"] = [
            (origin + numpy.dot(multi_indices, axes_rearranged)) * unit_conversion
            for multi_indices in numpy.indices(nrs_rearranged).reshape(3, -1).T
        ]
    f.close()
    return result


def read_charges_cube(
    filename,
    match_word_order=False,
    match_axis_order=True,
    add_nuclear_charges=False,
    force_angstroms=False,
    invert_charge_data=False,
    rescale_charges=True,
    total_charge=0,
    nr_return=None,
    density=False,
):
    """Read in a Gaussian-Cube file.
    
    Args:
        filename: (string) the name of the cube file
        match_word_order: (bool) if True, try to find out the order of inner,
            outer and middle to guess the correct volumetric data. The words
            outer, inner and middle have to be pressent in a certain order. For
            example the order OUTER X, INNER Y, MIDDLE Z will result in x, y, z being
            the outer, inner and middle loops, respectively. Per default,
            outer, middle and inner loop are x,y and z, respectively.
        match_axis_order: (bool) if True, try to find out the order of X, Y and
            Z coordinates of the volumetric data. The letters X, Y and Z have
            to be present in a certain order. Default: see match_word_order
        add_nuclear_charges: (bool) if True, nuclear charges and coordinates
            will be the first entries in the returned lists.
        force_angstroms: (bool) if this is True, enforce everything to be considered to
            be given in Angstroms
        invert_charge_data: (bool) if True, invert the volumetric charge data.
            Nuclear charges are always positive and volumetric data is taken as is.
        rescale_charges: (bool) the sum of atomic charges and the sum of
            nuclear charges have to match. If this is True and the charges
            don't match, rescale all volumetric data linearly so that they do.
            Only makes sense if add_nuclear_charges == True
        total_charge: (float) the total charge of the molecule to properly
            rescale the electronic charges
        nr_return: (list or None) if a variable of type list is given, append
            to it the number of atoms and the number of volumetric entries
        density: (bool) if True, return the density at the center of the voxel
            instead of the product of the density and the voxel's volume

    Returns:
        tuple of 2 lists. The first list contains 3-element lists of floats
        that are the Cartesian coordinates of the volumetric data values. The
        second list is a flat list of floats that contains the associated
        volumetric data. Nuclear charges can also be included.
    """
    f = open(filename, "r")
    next(f)
    line = next(f).rstrip()
    # per default, outer, middle and inner loop are x,y and z, respectively
    matching_wordperm = ["outer", "middle", "inner"]
    matching_xyzperm = ["x", "y", "z"]
    match_dict = {"outer": 0, "inner": 0, "middle": 0}
    xyz_dict = {"x": 0, "y": 1, "z": 2}
    # try to find out the order of the words outer, inner and middle as well as the order of x, y and z
    if match_word_order:
        if None in (
            re.search(regex, line, re.IGNORECASE)
            for regex in ["\\bouter\\b", "\\binner\\b", "\\bmiddle\\b"]
        ):
            raise ValueError(
                "You requested to match the word order against the second line but I cannot find the words outer, inner, middle there."
            )
        else:
            matching_wordperm = []
            for perm in itertools.permutations(["outer", "inner", "middle"]):
                if re.match(
                    ".*\\b" + "\\b.*\\b".join(perm) + "\\b.*", line, re.IGNORECASE
                ):
                    matching_wordperm = perm
                    break
    if match_axis_order:
        if None in (
            re.search(regex, line, re.IGNORECASE)
            for regex in ["\\bx\\b", "\\by\\b", "\\bz\\b"]
        ):
            raise ValueError(
                "You requested to match the order against the second line but I cannot find the words x, y, z there."
            )
        else:
            matching_xyzperm = []
            for perm in itertools.permutations(["x", "y", "z"]):
                if re.match(
                    ".*\\b" + "\\b.*\\b".join(perm) + "\\b.*", line, re.IGNORECASE
                ):
                    matching_xyzperm = perm
                    break
    for word, letter in zip(matching_wordperm, matching_xyzperm):
        match_dict[word] = xyz_dict[letter]
    # preallocate stuff
    axes = [None] * 3
    nrs = [None] * 3
    # prepare array for unit conversion
    # default in file is atomic units but everything will be transformed to Angstroms
    unit_conversion = numpy.array([0.5291772488] * 3, dtype=float)
    # read in the header
    # line with origin and number of atoms
    line = next(f).rstrip().split()
    nr_atoms = int(line[0])
    origin = numpy.array(list(map(float, line[1:4])))
    # the three lines with the voxel sizes and numbers of entries
    for c in range(3):
        line = next(f).rstrip().split()
        nrs[c] = abs(int(line[0]))
        if int(line[0]) < 0 or force_angstroms:
            unit_conversion[c] = 1.0
        axes[c] = list(map(float, line[1:4]))
    # the units of the origin have to be converted as well! doh...
    oring = origin * unit_conversion
    # convert to numpy arrays
    axes = numpy.array(axes)
    nrs = numpy.array(nrs)
    # this is the volume of one voxel
    # which will be used to convert charge density to charge
    # so it will seem as if there were a point charge at the center
    # of the voxel containing the whole charge inside that voxel
    if density:
        volume = 1.0
    else:
        volume = numpy.linalg.det(unit_conversion * axes.transpose())
    # read in volumetric data
    if add_nuclear_charges:
        sum_nuclear_charges = total_charge
        charges = numpy.zeros((nrs[0] * nrs[1] * nrs[2] + nr_atoms), dtype=float)
        coordinates = numpy.zeros((nrs[0] * nrs[1] * nrs[2] + nr_atoms, 3), dtype=float)
        for count in range(nr_atoms):
            line = next(f).rstrip().split()
            # nuclear charges have to have the opposite sign as electronic charges
            charges[count] = int(line[0])
            sum_nuclear_charges += charges[count]
            # the first coloumn contains the atomic charge and the second is undefined
            # so the last 3 contain the information I need`
            coordinates[count] = list(map(float, line[2:5]))
        count = nr_atoms
    else:
        # skip lines of atomic positions
        charges = numpy.zeros((nrs[0] * nrs[1] * nrs[2]), dtype=float)
        coordinates = numpy.zeros((nrs[0] * nrs[1] * nrs[2], 3), dtype=float)
        count = 0
        for i in range(nr_atoms):
            next(f)
    sum_electronic_charges = 0
    for l in f:
        for e in map(float, l.rstrip().split()):
            charges[count] = e * volume
            sum_electronic_charges += charges[count]
            count += 1
    if add_nuclear_charges and (rescale_charges or invert_charge_data):
        is_nucleus = numpy.zeros(charges.shape, dtype=bool)
        is_nucleus[:nr_atoms] = numpy.ones((nr_atoms), dtype=bool)
    if add_nuclear_charges and rescale_charges:
        electronic_charge_rescale_factor = sum_nuclear_charges / sum_electronic_charges
        charges = (
            charges * is_nucleus
            + charges * numpy.logical_not(is_nucleus) * electronic_charge_rescale_factor
        )
    if invert_charge_data:
        if add_nuclear_charges:
            charges = charges * is_nucleus + charges * numpy.logical_not(is_nucleus) * (
                -1
            )
        else:
            charges *= -1
    axes_rearranged = numpy.zeros(axes.shape, dtype=float)
    axes_rearranged[0] = axes[match_dict["outer"]]
    axes_rearranged[1] = axes[match_dict["middle"]]
    axes_rearranged[2] = axes[match_dict["inner"]]
    nrs_rearranged = numpy.zeros(nrs.shape, dtype=int)
    nrs_rearranged[0] = nrs[match_dict["outer"]]
    nrs_rearranged[1] = nrs[match_dict["middle"]]
    nrs_rearranged[2] = nrs[match_dict["inner"]]
    unit_conversion_rearranged = numpy.zeros(unit_conversion.shape, dtype=float)
    unit_conversion_rearranged[0] = unit_conversion[match_dict["outer"]]
    unit_conversion_rearranged[1] = unit_conversion[match_dict["middle"]]
    unit_conversion_rearranged[2] = unit_conversion[match_dict["inner"]]
    if add_nuclear_charges:
        count = nr_atoms
    else:
        count = 0
    # o, m, i stand for outer, inner and middle, respectively
    # All the 4 variants do the exact same thing, but the last is approx. twice as fast as the first
    # Variant 1
    # for multi_indices in (numpy.array((o,m,i)) for o in range(nrs_rearranged[0]) for m in range(nrs_rearranged[1]) for i in range(nrs_rearranged[2])):
    #    position=(origin+numpy.dot(multi_indices,axes_rearranged))*unit_conversion_rearranged
    #    coordinates[count]=position
    #    count+=1
    # Variant 2
    # for o in range(nrs_rearranged[0]):
    #    for m in range(nrs_rearranged[1]):
    #        for i in range(nrs_rearranged[2]):
    #            position=(origin+numpy.dot(numpy.array((o,m,i)),axes_rearranged))*unit_conversion_rearranged
    #            coordinates[count]=position
    #            count+=1
    # Variant 3
    # for multi_indices in numpy.indices(nrs_rearranged).reshape(3,-1).T:
    #    position=(origin+numpy.dot(multi_indices,axes_rearranged))*unit_conversion_rearranged
    #    coordinates[count]=position
    #    count+=1
    # variant 4
    coordinates[count:] = [
        (origin + numpy.dot(multi_indices, axes_rearranged))
        * unit_conversion_rearranged
        for multi_indices in numpy.indices(nrs_rearranged).reshape(3, -1).T
    ]
    coordinates[:count] = coordinates[:count] * unit_conversion
    close(f)
    if not nr_return == None:
        nr_return.append(nr_atoms)
        nr_return.append(len(charges) - nr_atoms)
    return coordinates, charges


def _string_to_boolean(string):
    """Convert a string to boolean, case-insensitively.

    Args:
        string: (string) string representation of a boolean value

    Raises:
        TypeError.
    """
    v = str(string).lower()
    if v in ["true", "false"]:
        return v == "true"
    else:
        raise TypeError("Not a boolean: %s" % (string))


# taken from http://stackoverflow.com/questions/2885190/using-pythons-configparser-to-read-a-file-without-section-name
# and modified to be less elaborate
class SectionlessConfigParser(ConfigParser.ConfigParser):
    """Extends ConfigParser to allow files without sections.

    This is done by wrapping read files and prepending them with a placeholder
    section, which defaults to '__DEFAULT__'.

    Create an object of this class using the function
    ManipulateAggregates.collection.read.read_config_file
    """

    def __init__(self, nocase=False, sep=None, *args, **kwargs):
        """Constructor.
        
        No arguments required. Do not use directly. Use
        ManipulateAggregates.collection.read.read_config_file instead
        """
        self.nocase = nocase
        self.sep = sep
        ConfigParser.ConfigParser.__init__(self, *args, **kwargs)

    def _readfp(self, fp, *args, **kwargs):
        """Open the config file and read it in.

        Args:
            fp: (string) config file name

        Returns:
            object of class ConfigParser.ConfigParser
        """
        if self.nocase:
            lower1st = lambda l: (
                i.lower() if c == 0 else i for i, c in list(zip(l, range(len(l))))
            )
            if self.sep is not None:
                sep = self.sep
            else:
                sep = "="
            translate = lambda l: "=".join(lower1st(l.split(sep, 1)))
        else:
            if self.sep is not None:
                translate = lambda l: "=".join(l.split(self.sep, 1))
            else:
                translate = lambda l: l
        with open(fp, "r") as stream:
            lines = (translate(l) for l in stream.readlines())
            fakefile = StringIO(tounicode("[__DEFAULT__]\n" + "\n".join(lines)))
        return readfp(self, fakefile, *args, **kwargs)

    def _convert(self, func, name, *args, **kwargs):
        """Convert an entry using a function.

        Args:
            func: (function) this function is applied to the value associated
                with the key name
            name: (string) key whose associated value will be returned

        Returns:
            the converted value

        Raises:
            TypeError.
        """
        try:
            v = self.get("__DEFAULT__", name, *args, **kwargs)
        except ConfigParser.InterpolationMissingOptionError as e:
            errorstring = (
                "Keyword '%s' requested (which defaults to the value of keyword '%s') but no value could be found."
                % (e[0], e[3])
            )
            raise NoOptionInConfigFileError(errorstring)
        except ConfigParser.NoOptionError as e:
            errorstring = "Keyword '%s' requested but no value could be found." % (e[0])
            raise NoOptionInConfigFileError(errorstring)
        try:
            return func(v)
        except TypeError as e:
            raise TypeError(
                "Value associated with keyword '%s' is of wrong type." % name, e
            )

    def get_int(self, name, *args, **kwargs):
        """Get an integer.

        Returns:
            the value associated with the key name converted to integer.

        Raises:
            ValueError.
        """
        try:
            return self._convert(int, name, *args, **kwargs)
        except ValueError as e:
            raise ValueError(
                "Value associated with keyword '%s' could not be converted to int."
                % name,
                e,
            )

    def get_float(self, name, *args, **kwargs):
        """Get a floating point value.

        Returns:
            the value associated with the key name converted to float.

        Raises:
            ValueError.
        """
        try:
            return self._convert(float, name, *args, **kwargs)
        except ValueError as e:
            raise ValueError(
                "Value associated with keyword '%s' could not be converted to float."
                % name,
                e,
            )

    def get_boolean(self, name, *args, **kwargs):
        """Get a boolean.

        Returns:
            the value associated with the key name converted to bool.

        Raises:
            ValueError.
        """
        try:
            return self._convert(_string_to_boolean, name, *args, **kwargs)
        except ValueError as e:
            raise ValueError(
                "Value associated with keyword '%s' could not be converted to boolean."
                % name,
                e,
            )

    def _allitems(self, *args, **kwargs):
        """Get all items in the config file.

        Returns:
            a list of config options (strings) in the config file
        """
        try:
            return self.items("__DEFAULT__", *args, **kwargs)
        except ConfigParser.InterpolationMissingOptionError as e:
            errorstring = (
                "Keyword '%s' requested (which defaults to the value of keyword '%s') but no value could be found."
                % (e[0], e[3])
            )
            raise NoOptionInConfigFileError(errorstring)

    def get_str(self, name, *args, **kwargs):
        """Get a string.

        Returns:
            the value associated with the key name converted to string.
        """
        return self._convert(str, name, *args, **kwargs)

    def check_against(self, options):
        """Check options against the content of the config file.
        
        Args:
            options: (list of strings) options that are expected in the config
                file. options may contain more entries that what are
                expected but any option in the config file but not in options
                will be returned.

        Returns:
            a list of strings that contain lines with keyword and value of
            those keywords that are not in options but in the config file
        """
        return ["%-20s = %s" % o for o in self._allitems() if not o[0] in options]


def read_config_file(filename, defaults=None, nocase=False, sep=None):
    """Read in a section-less config file.

    Use methods "get_str('name')" of the returned object to get a string object
    corresponding to the key "name".  Other functions defined: get_int,
    get_float, get_boolean, _allitems (return all items as strings in a
    dictionary). Types appropriate to the name will be returned.

    Args:
        filename: (string) the name of the config file to read in
        defaults: (dictionary) a dictionary providing default values
        nocase: (bool) whether or not to ignore the case of the keywords
        sep: (string) if using a non-standard cfg-file separator, specify it here
    
    Returns:
        an object of ManipulateAggregates.collection.read.SectionlessConfigParser
    """
    # info on how to do this has been taken from:
    # http://stackoverflow.com/questions/2885190/using-pythons-configparser-to-read-a-file-without-section-name
    parser = SectionlessConfigParser(defaults=defaults, nocase=nocase, sep=sep)
    parser._readfp(filename)
    return parser


def _gen_triples_off(iterable, check_first=False, convert_func=lambda x: x):
    """Generate tuples of 3 elements from an iterable.

    If the number of elements in iterable is not divisible by 3, skip the
    last up to 2 elements.

    Args:
        iterable: (iterable) flat iterable
        check_first: (bool) if True, assume that each tuple is preceeded by
            the number of elements. An exception is raised if that number is
            not equal to 3.
        convert_func: (function) apply this function to each element of the
            iterable prior to returning the triples.

    Returns:
        a generator that yiels 3-element tuples

    Raises:
        ValueError.
    """
    it = (i for i in iterable)
    try:
        while True:
            if check_first:
                chk = int(next(it))
                if chk != 3:
                    print(chk)
                    raise ValueError("Face with 4 vertices detected.")
            a = convert_func(next(it))
            b = convert_func(next(it))
            c = convert_func(next(it))
            yield (a, b, c)
    except StopIteration:
        return


def read_off(filename):
    """Read an OFF file.

    Args:
        filename: (string) the name of the OFF file

    Returns:
        a list of 3-element tuples. Each of these tuples contains the Cartesian
            coordinates of a vertex. Each of the lists describes a face of the
            surface that is read in.
    """
    f = open(filename, "r")
    # read in data
    lines = [l.rstrip().split() for l in f]
    # flatten list since linebreaks don't matter in that filetype
    lines = [e for l in lines for e in l]
    close(f)
    if lines[0].lower() != "off":
        raise ValueError("The first entry has to be OFF.")
    nr_vertices = int(lines[1])
    nr_faces = int(lines[2])
    nr_half = int(lines[3])
    if nr_half != 0:
        raise ValueError("Number of half-edges is unequal zero. This is not supported.")
    vertices = list(
        _gen_triples_off(lines[4 : 4 + 3 * nr_vertices], convert_func=float)
    )
    face_indices = _gen_triples_off(lines[4 + 3 * nr_vertices :], True, int)
    faces = [tuple(vertices[i] for i in face) for face in face_indices]
    return faces
