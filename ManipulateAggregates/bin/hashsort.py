"""Executable to sort files into hashed subdirectories.

See the variable hashsort.HELPTEXT for more details.
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

# import sys,os,logging
# logfile = os.getenv("MALOGFILE",None)
# loglevel = getattr(logging,os.getenv("MALOGLEVEL","INFO").upper())
# logging.basicConfig(filename=logfile,level=loglevel)
# logger = logging.getLogger("hashsort")
#
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
from ManipulateAggregates.collection import hashIO

# help message text
HELPTEXT = """This script allows to copy files into subdirectories depending on
the hashes of their names. Usage:

hashsort.py DIRECTORY [REGEX]

The mandatory first argument DIRECTORY must be a directory that already exists.
All files in this directory whose names match the regular expression provided
as the second, optional argument (default: '^[1-9][0-9]*_out.dx$') will be sorted
into subdirectories. This will use the MD5 hasing algorithm.
"""

# default process name
PROCNAME = "HashSort"
try:
    from FireDeamon import set_procname

    set_procname(PROCNAME)
except ImportError:
    set_procname = lambda s: None

# default regular expression
REGEX = "^[1-9][0-9]*_out.dx$"


def entrypoint():
    for arg in sys.argv:
        if arg == "--help" or arg == "-h":
            print(HELPTEXT)
            exit(0)
    if len(sys.argv) < 3:
        raise ValueError("Not enough arguments given.")
    dir = sys.argv[1]
    if not os.path.isdir(dir):
        raise IOError("First argument must be a directory")
    dxregex = re.compile(REGEX)
    if len(sys.argv) == 3:
        dxregex = re.compile(sys.argv[2])
    count = 0
    for f in os.listdir(dir):
        if os.path.isfile(dir + os.sep + f) and dxregex.match(f) is not None:
            count += 1
            print("%s -> %s" % (dir + os.sep + f, hashIO.hashpath(dir + os.sep + f)))
            hashIO.exists(dir + os.sep + f, nulldepth=True)
    print("Sorted %d files into hashed directories." % (count))


if __name__ == "__main__":
    entrypoint()
