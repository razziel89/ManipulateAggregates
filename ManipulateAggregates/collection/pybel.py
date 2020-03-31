"""
pybel - A Cinfony module for accessing Open Babel

Global variables:
  ob - the underlying SWIG bindings for Open Babel
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  forcefields - a list of supported forcefields

Modified to be a part of ManipulateAggregates
"""

# Copyright (c) 2008-2012, Noel O'Boyle; 2012, Adria Cereto-Massague
# Copyright (C) 2020 by Torsten Sachse
# All rights reserved.
#
# This file is part of ManipulateAggregates.
# The contents of this file are covered by the terms of the GPL v2 license,
# see the file LICENSE_GPLv2 in the root if the repository for the license.
#
# ManipulateAggregates is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version. This file, however, is covered by the
# the terms of the GNU General Public License as published by
# the Free Software Foundation in version 2 of the License.
#
# ManipulateAggregates is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ManipulateAggregates. If not, see <http://www.gnu.org/licenses/>.


import sys
import os.path

import openbabel as ob

_obfuncs = _obconsts = ob


def _formatstodict(list):
    broken = [
        x.replace("[Read-only]", "").replace("[Write-only]", "").split(" -- ")
        for x in list
    ]
    broken = [(x, y.strip()) for x, y in broken]
    return dict(broken)


def _getplugins(findplugin, names):
    return dict([(x, findplugin(x)) for x in names if findplugin(x)])


def _getpluginnames(ptype):
    plugins = ob.vectorString()
    ob.OBPlugin.ListAsVector(ptype, None, plugins)
    return [x.split()[0] for x in plugins]


_obconv = ob.OBConversion()
_builder = ob.OBBuilder()

# A dictionary of supported input formats
informats = _formatstodict(_obconv.GetSupportedInputFormat())
# A dictionary of supported output formats
outformats = _formatstodict(_obconv.GetSupportedOutputFormat())

# A list of supported forcefields
forcefields = [_x.lower() for _x in _getpluginnames("forcefields")]
_forcefields = _getplugins(ob.OBForceField.FindType, forcefields)


def readfile(format, filename, opt=None):
    """Iterate over the molecules in a file.

    Args:
        format: see the informats variable for a list of available input formats
        filename: tha name of the file from which to read the molecular data
        opt: a dictionary of format-specific options For format options with no
            parameters, specify the value as None.

    You can access the first molecule in a file using the next() function:
    >>> mol = next(readfile("smi", "myfile.smi"))

    You can make a list of the molecules in a file using:
    >>> mols = list(readfile("smi", "myfile.smi"))

    You can iterate over the molecules in a file as shown in the following code snippet:
    >>> atomtotal = 0
    >>> for mol in readfile("sdf", "head.sdf"):
    >>>     atomtotal += len(mol.atoms)
    """
    if opt is None:
        opt = {}
    obconversion = ob.OBConversion()
    formatok = obconversion.SetInFormat(format)
    for k, v in opt.items():
        if v is None:
            obconversion.AddOption(k, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(k, obconversion.INOPTIONS, str(v))
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" % format)
    if not os.path.isfile(filename):
        raise IOError("No such file: '%s'" % filename)

    def filereader():
        obmol = ob.OBMol()
        notatend = obconversion.ReadFile(obmol, filename)
        while notatend:
            yield Molecule(obmol)
            obmol = ob.OBMol()
            notatend = obconversion.Read(obmol)

    return filereader()


class Outputfile(object):
    """Represent a file to which *output* is to be sent.

    Although it's possible to write a single molecule to a file by calling the write()
    method of a molecule, if multiple molecules are to be written to the same file you
    should use the Outputfile class.

    Args:
        format: see the outformats variable for a list of available output formats
        filename
        overwrite: if the output file already exists, should it be overwritten? (default
            is False)
        opt: a dictionary of format-specific options For format options with no
            parameters, specify the value as None.

    Methods:
       write(molecule)
       close()
    """

    def __init__(self, format, filename, overwrite=False, opt=None):
        if opt is None:
            opt = {}
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError(
                "%s already exists. Use 'overwrite=True' to overwrite it."
                % self.filename
            )

        self.obConversion = ob.OBConversion()
        formatok = self.obConversion.SetOutFormat(self.format)
        if not formatok:
            raise ValueError("%s is not a recognised Open Babel format" % format)
        if filename and filename.split(".")[-1] == "gz":
            self.obConversion.AddOption("z", self.obConversion.GENOPTIONS)
        for k, v in opt.items():
            if v is None:
                self.obConversion.AddOption(k, self.obConversion.OUTOPTIONS)
            else:
                self.obConversion.AddOption(k, self.obConversion.OUTOPTIONS, str(v))
        self.total = 0  # The total number of molecules written to the file

    def write(self, molecule):
        """Write a molecule to the output file.

        Args:
            molecule
        """
        if not self.filename:
            raise IOError("Outputfile instance is closed.")

        if self.total == 0:
            self.obConversion.WriteFile(molecule.OBMol, self.filename)
        else:
            self.obConversion.Write(molecule.OBMol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.obConversion.CloseOutFile()
        self.filename = None


class Molecule(object):
    """Represent a Pybel Molecule.

    Args:
        OBMol: an Open Babel OBMol

    Methods:
       addh(), localopt(), make3D(), removeh(), write()

    The underlying Open Babel molecule can be accessed using the attribute:
       OBMol
    """

    def __init__(self, OBMol):
        self.OBMol = OBMol

    @property
    def charge(self):
        return self.OBMol.GetTotalCharge()

    @property
    def conformers(self):
        return self.OBMol.GetConformers()

    @property
    def dim(self):
        return self.OBMol.GetDimension()

    @property
    def energy(self):
        return self.OBMol.GetEnergy()

    @property
    def exactmass(self):
        return self.OBMol.GetExactMass()

    @property
    def formula(self):
        return self.OBMol.GetFormula()

    @property
    def molwt(self):
        return self.OBMol.GetMolWt()

    @property
    def spin(self):
        return self.OBMol.GetTotalSpinMultiplicity()

    def _gettitle(self):
        return self.OBMol.GetTitle()

    def _settitle(self, val):
        self.OBMol.SetTitle(val)

    title = property(_gettitle, _settitle)

    def write(self, format="smi", filename=None, overwrite=False, opt=None):
        """Write the molecule to a file or return a string.

        If a filename is specified, the result is written to a file. Otherwise, a string
        is returned containing the result.

        To write multiple molecules to the same file you should use the Outputfile
        class.

        Args:
            format: see the informats variable for a list of available output formats
                (default is "smi")
            filename: default is None
            overwite: if the output file already exists, should it be overwritten?
                (default is False)
            opt: a dictionary of format specific options For format options with no
                parameters, specify the value as None.

        """
        if opt is None:
            opt = {}
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError("%s is not a recognised Open Babel format" % format)
        if filename and filename.split(".")[-1] == "gz":
            obconversion.AddOption("z", self.obConversion.GENOPTIONS)
        for k, v in opt.items():
            if v is None:
                obconversion.AddOption(k, obconversion.OUTOPTIONS)
            else:
                obconversion.AddOption(k, obconversion.OUTOPTIONS, str(v))

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError(
                    ("%s already exists. Use 'overwrite=True' to " "overwrite it.")
                    % filename
                )
            obconversion.WriteFile(self.OBMol, filename)
            obconversion.CloseOutFile()
        else:
            return obconversion.WriteString(self.OBMol)

    def localopt(self, forcefield="mmff94", steps=500):
        """Locally optimize the coordinates.

        If the molecule does not have any coordinates, make3D() is called before the
        optimization. Note that the molecule needs to have explicit hydrogens. If not,
        call addh().

        Args:
            forcefield: default is "mmff94". See the forcefields variable for a list of
                available forcefields.
            steps: default is 500
        """
        forcefield = forcefield.lower()
        if self.dim != 3:
            self.make3D(forcefield)
        ff = _forcefields[forcefield]
        success = ff.Setup(self.OBMol)
        if not success:
            raise RuntimeError("Local optimization with forcefield unsuccessful")
        ff.SteepestDescent(steps)
        ff.GetCoordinates(self.OBMol)

    def make3D(self, forcefield="mmff94", steps=50):
        """Generate 3D coordinates.

        Once coordinates are generated, hydrogens are added and a quick local
        optimization is carried out with 50 steps and the MMFF94 forcefield. Call
        localopt() if you want to improve the coordinates further.

        Args:
            forcefield: default is "mmff94". See the forcefields variable for a list of
                available forcefields.
            steps: default is 50
        """
        forcefield = forcefield.lower()
        _builder.Build(self.OBMol)
        self.addh()
        self.localopt(forcefield, steps)

    def addh(self):
        """Add hydrogens."""
        self.OBMol.AddHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.OBMol.DeleteHydrogens()

    def __str__(self):
        return self.write()
