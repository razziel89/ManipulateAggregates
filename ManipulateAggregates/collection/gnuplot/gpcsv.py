"""Convert data from CSV files to a format that the package gnuplot can handle.

CSV files can contain data with a lot of different delimiters. This file
provides functionality to test a variety of combinations of delimiters and
choses the first paring that succeeds.

Attributes:
    delimiters (list of tuples of characters): the set of delimites that will
        be checked against. [(",","."),(";",","),(" ","."),("\t","."),(" ",",")]
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

import csv
import re
import io
import copy

from . import postprocess as pp

delimiters = [(",", "."), (";", ","), (" ", "."), ("\t", "."), (" ", ",")]


def gpcsv(
    filename,
    GP,
    xcol=0,
    ycols=None,
    delimiter=None,
    separator=None,
    default=None,
    ppargs=None,
    args=None,
):
    """Convert the data in a CSV-file to the Gnuplot format.

    A certain set of delimiters is checked and the first combination that
    yields valid floats is chosen (if no delimiters are provided).  Delimiter
    and separator default to a valid combination from the following list:
    [(",","."),(";",","),(" ","."),("\\t","."),(" ",",")]

    Args:
        filename: (string) The name of the CSV-file to be converted.
        GP: (gnuplot) The object of the gnuplot class that should plot the data.
        xcol: (int) The coloumn (starting at 0) in the CSV file that contains
            the data to be used on the x-axis. Optional, default: 0
        ycols: (list of int) The coloumns in the CSV file that contain the data
            sets to be used on the y-axis. Optional, default: None
        delimiter: (string) The coloumn delimiter for the CSV-file.
            Optional, default: try all from default combinations.
        separator: (string) The decimal separator for the CSV-file.
            Optional, default: try all from default combinations.
        default: (dict) Copy all configuration parameters from this dictionary
            (using copy.copy) to get default values. Optional, default: empty
        ppargs: (list of strings) If given, apply functions known by that
            name in order to the extracted data. Optional, default: None
        args: (list of lists of arguments) If given, pass these
            arguments in order to the functions in postprocess.
            Optional, default: None
    """
    if delimiter is not None and separator is not None:
        delim, sep = delimiter, separator
    else:
        delim, sep = (None, None)
        for d, s in delimiters:
            with io.open(filename, "r", newline="") as csvfile:
                spamreader = csv.reader(csvfile, delimiter=d)
                for row in spamreader:
                    newrow = [r.replace(s, ".") for r in row if len(r) > 0]
                    if len(newrow) > 0:
                        try:
                            list(map(float, newrow))
                        except ValueError:
                            # if the above map command did not work on a line that
                            # contained only numbers and the delimiters, this pattern
                            # does not apply
                            if re.match("[0-9" + s + d + "]", d.join(newrow)):
                                break
                # the else branch is executed if the end of the for loop has been
                # reached break avoids this to be reached and, hence, else branch is not
                # executed
                else:
                    delim, sep = (d, s)
                    break
        if delim is None or sep is None:
            raise IOError(
                "Could not correctly parse the given CSV-file %s with any of the available delimiters."
                % (filename)
            )
    data = []
    nr_fields = float("inf")
    with io.open(filename, "r", newline="") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=d)
        for row in spamreader:
            newrow = [r.replace(s, ".") for r in row if len(r) > 0]
            if len(newrow) > 0:
                try:
                    l = list(map(float, newrow))
                    if len(l) < nr_fields:
                        nr_fields = len(l)
                    data.append(l)
                except ValueError:
                    # if the above map command did not work on a line that contained
                    # only numbers and the delimiters, this pattern does not apply
                    if re.match("[0-9" + s + d + "]", d.join(newrow)):
                        raise IOError(
                            "Could not correctly parse the given CSV-file %s with the given delimiters %s,%s."
                            % (filename, delim, sep)
                        )
    if ycols is None:
        ycols = list(range(1, nr_fields))
    olddata = copy.deepcopy(data)
    if ppargs is not None and args is not None:
        data = pp.apply(data, ppargs, args, xcol, ycols)
    tempfile = GP.data_to_file(data, formatstring=None, delete=True)
    if default is not None and isinstance(default, dict):
        dictionary = copy.copy(default)
    else:
        dictionary = {}
    if not "title" in dictionary:
        dictionary["title"] = filename.replace("_", " ")
    dictionary["type"] = "file"
    dictionary["file"] = tempfile
    return dictionary
