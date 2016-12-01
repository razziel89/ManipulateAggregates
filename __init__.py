"""The top-level module ManipulateAggregates.

Please refer to the submodules' documentations for more details.

Submodules:
  - collection:
      read and write several file formats used in computational chemistry,
      control OpenGL from Python energyscan: estimate energetically
      favourable aggregate (dimers and higher ones) geometries in a
      three-step procedure
  - aggregate:
      manipulate (internal) degrees of freedom of molecules and aggregates
      and visualize distributions of electrostatic potentials and electron
      densities
  - orbitalcharacter:
      compute electrostatic potentials and electron densities from quantum
      chemical orbitals

@package ManipulateAggregates
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

if __name__ == "__main__":
    raise Exception("Do not run ManipulateAggregates directly.")

import sys,os,logging
logfile = os.getenv("MALOGFILE",None)
loglevel = getattr(logging,os.getenv("MALOGLEVEL","WARNING").upper())
logging.basicConfig(filename=logfile,level=loglevel)
logger = logging.getLogger(__name__)
if logfile is None:
    logger.info("Set the environment variable MALOGFILE to log errors to a file of that name.")

from . import collection as collection
from . import energyscan as energ
from . import aggregate as aggregate
from . import orbitalcharacter as orbitalcharacter
