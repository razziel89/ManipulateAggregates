# -*- coding: utf-8 -*-
"""Create low-energy aggregates from molecules.

This is the submodule energyscan. It estimates energetically favourable
aggregate (dimers and higher ones) geometries in a three-step procedure.

1. create a huge number of aggregates:
   - implemented in ManipulateAggregates.energyscan.scan
   - Keep a central molecule fixed and
   - move another molecule around the central one (varying orientations) and
   - evaluate energy for every aggregate created (so far: using force fields).
   - molecules can be replaced by entire aggregates
2. search for local energy minima among the aggregates created in the previous step
   - implemented in ManipulateAggregates.energyscan.minimasearch
3. screen all of those local energy minimum structures to obtain a highly diverse set
   - implemented in ManipulateAggregates.energyscan.similarityscreening

It requires a custom and reduced version of OpenBabel
(https://github.com/razziel89/MaAg-bel) including its Python bindings. It also requires
the Python package FireDeamon (https://github.com/razziel89/libfiredeamon).

It currently uses OpenBabel's force fields to estimate aggregate energies.
Please see the global help text variable energyscan.LONGHELPTEXT for a more
detailed description.

Parallelization is supported for some subsubmodules. Set the environment
variable OMP_NUM_THREADS to the value you want to use. Only single-node
parallelization is supported as of now.

Details:

The generation of low-energy aggregates is a three-step procedure. The first
step is a volumetric energy scan: first, two regular three-dimensional grids
are generated, the molecules’ main axes are aligned with the three Cartesian
axes and they are centered at the origin. Then, one of the molecules is kept
fixed while the other’s center is moved to every point on the first, i.e.,
spatial grid and the aggregate’s energy is evaluated. Afterwards, the second
molecule is rotated according to Euler angles around its main axes defined by
the second, i.e., angular grid and the scan is repeated for every such
orientation. I greatly improved the speed of this brute-force approach by
excessive screening of vdW-overlaps and aggregates, whose vdW surfaces are
farther apart than a given distance (here: 2.5 Å). I remark that this screening
has to be performed for every angular arrangement independently and only very
few points can be screened right at the start, i.e., those inside and slightly
outside the vdW-surface of the central molecule.

The second step is the search for geometries corresponding to local energy
minima for each angular arrangement: I define a local minimum spatial point as
such a point, whose associated aggregate has an energy lower than that of every
other point in its neighborhood. Given a definition for “neighborhood”, a
neighbor list is generated associating each point with its neighbors allowing
for treating irregular grids just like regular ones. This list is then used to
determine all local energy minima by simple comparison of the stored volumetric
data.

The third step is the selection of the desired number of highly diverse
low-energy aggregate geometries: first, the programme discards all those
remaining geometries whose energy, possibly obtained using a different method
than before, is higher than a given cutoff. Geometrical similarity is then
estimated by computing the RMSD between two aggregates after they were again
aligned to the origin and the three Cartesian axes. In order to minimize the
required number of RMSD computations, I use a tree-based algorithm as
implemented in OpenBabel. This algorithm, however, is designed to quickly
extract such geometries from a given set, whose RMSDs with respect to each
other are all larger than a given cutoff, without storing pairwise RMSDs. In
order to obtain a small number of highly diverse dimers, the aforementioned
cutoff is increased step by step, thereby successively reducing the number of
dimers until there is at most the desired number of them left.

The aggregates’ energies were evaluated using force-fields as implemented in
the OpenBabel chemical toolbox. Technically, every energy evaluation method
could be used for this algorithm. However, the great advantage of force-fields
is the high speed of energy evaluations as compared to, for example, density
functional theory calculations. The great disadvantage of force-fields is
undoubtedly that some bonds might not be parameterized correctly for newly
synthesized molecules. Since this algorithm evaluates differences only in
intermolecular interactions, I believe it to be accurate enough to determine
starting structures for further treatment as the here described algorithm is
only the first step in the computational generation of higher aggregates
(described elsewhere, not yet published). In order to avoid the generation of
false low-energy local minimum aggregates, the resulting structures should be
quantum chemically optimized in a final step. 
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
from . import ansilliary
from . import minimasearch
from . import scan
from . import similarityscreening
