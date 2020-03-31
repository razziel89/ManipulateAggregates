.. ManipulateAggregates documentation master file, created by
   sphinx-quickstart on Thu Mar 26 12:47:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ManipulateAggregates's documentation!
================================================

.. toctree::
    :maxdepth: 3
    :caption: Contents:


Introduction
============

What you are currently viewing contains the documentation for
ManipulateAggregates, a set of Python scripts written to perform some tasks
related to what I did during my time as a PhD student that will be detailed in
this documentation. For any license-related information, please see the file
called COPYING and the header of each individual `*.py` file.

The module ManipulateAggregates consists of four submodules each of which
resides in its own directory. The following is a list of all four submodules (in
alphabetical order) including a synopsis of the functionality they provide.

1. collection:
    - read and write several file formats used in computational chemistry:
        - geometry: cube, molden, xyz
        - orbital data: molden
        - volumetric data: cube, dx, xyz
        - frequencies: aims, terachem
    - control OpenGL from Python:
        - draw (coloured) spheres and trimeshes
        - export to pov-file to render with PoVRay
    - prefix file names by auto-generated hashes to limit the number of files
      per directory (huge speed-ups for some file systems)
    - read section-less Python config files (with defaults and type checks)
    - most file types can also be read and written when gzipped
    - control gnuplot from Python

2. energyscan:
    - estimate energetically favourable aggregate (dimers and higher ones)
      geometries in a three-step procedure:

        1. create a huge number of aggregates:
            - Keep a central molecule fixed and
            - move another molecule around the central one (varying
              orientations) and
            - evaluate energy for every aggregate created (so far: using force
              fields).
            - molecules can be replaced by entire aggregates
        2. search for local energy minima among the aggregates created in the
           previous step
        3. screen all of those local energy minimum structures to obtain a
           highly diverse set

    - requires a slightly modified version of `OpenBabel_` including its Python
      bindings (creation of a differently named fork containing only the
      functionality required here is a work-in-progress, allows using both
      versions of OpenBabel simultaneously)
    - requires `FireDeamon_`

3. aggregate:
    - manipulate (internal) degrees of freedom of molecules and aggregates
    - (compute and) visualize distributions of electrostatic potentials and
      electron densities

        - empirical methods supported via `OpenBabel_`
        - methods based on results from ab-initio computations supported via
          `libFireDeamon_`
        - live visualization using OpenGL 
        - high-quality visualization using PoVRay
        - support for volumetric distributions
        - support for distributions on surfaces:

            - isosurfaces through volumetric data
            - (scaled) van-der-Waals surfaces around molecules
            - arbitrarily high degrees of discretization supported

        - computations require the submodule orbitalcharacter
        - example image for electrostatic potential below

    - estimate a molecule's HLB value
    - support for all file types supported by OpenBabel

4. orbitalcharacter:
    - compute electrostatic potentials and electron densities from quantum
      chemical orbitals

        - corrections for the limited precision of the input data can be applied
        - computations can be restricted to use only some of the available
          orbitals

    - estimate the character of an orbital by comparison with orbitals of known
      character
    - much of the functionality provided within the submodule orbitalcharacter
      can be sped up when using the C++-library `libFireDeamon_`

Example image for the electrostatic potential on a van-der-Waals surface:

.. figure:: /images/example_vdw_pot.png

.. _OpenBabel: https://github.com/razziel89/openbabel
.. _libFireDeamon: https://github.com/razziel89/libfiredeamon


Prerequisites
=============

You need the following dependencies to use this package:
    * GNU make
    * Python:
        * `Anaconda`_-based virtual environments highly recommended
        * including the pip package manager
        * including the packages setuptools and wheel

If you want all the functionality specified above, you also need:
    * `OpenBabel_` (slightly modified version compared to upstream)
    * `libFireDeamon_`
    * NumPy

.. _Anaconda: https://www.anaconda.com/
.. _OpenBabel: https://github.com/razziel89/openbabel
.. _libFireDeamon: https://github.com/razziel89/libfiredeamon

Installation
============

If you have all the aforementioned dependencies installed and you environment
activated, simply run `pip install ManipulateAggregates`.


The Python Package
==================

.. automodule:: ManipulateAggregates
    :members:
    :undoc-members:

Aggregate
---------

.. automodule:: ManipulateAggregates.aggregate
    :members:
    :undoc-members:

Collection
----------

.. automodule:: ManipulateAggregates.collection
    :members:
    :undoc-members:

Collection.gnuplot
------------------

.. automodule:: ManipulateAggregates.collection.gnuplot
    :members:
    :undoc-members:

Collection.opengl
-----------------

.. automodule:: ManipulateAggregates.collection.opengl
    :members:
    :undoc-members:

Collection.pybel
----------------

.. automodule:: ManipulateAggregates.collection.pybel
    :members:
    :undoc-members:

EnergyScan
----------

.. automodule:: ManipulateAggregates.energyscan
    :members:
    :undoc-members:

Orbitalcharacter
----------------

.. automodule:: ManipulateAggregates.orbitalcharacter
    :members:
    :undoc-members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Misc
====

- Author: Torsten Sachse
- Date: 2015-2020
- Version: 0.1
- License:
    GNU General Public License version 3, apart from
    `ManipulateAggregates/collection/pybel.py`, which is under the GNU General Public
    License version 2

