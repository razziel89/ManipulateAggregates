"""This sobmodule is an aggregator for auxilliary functionality.

It provides functionality to:
  - read and write several file formats used in computational chemistry
  - control OpenGL from Python
  - prefix file names by auto-generated hashes to limit the number of files
    per directory (huge speed-ups for some file systems)
  - control gnuplot from Python

@package ManipulateAggregates.collection
"""

#This file is part of ManipulateAggregates.
#
#Copyright (C) 2016 by Torsten Sachse
#
#ManipulateAggregates is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ManipulateAggregates is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.

from . import gnuplot as gnuplot
from . import hashIO  as hashIO
from . import opengl  as opengl
from . import read    as read
from . import write   as write
from . import p2p3IO  as p2p3IO
