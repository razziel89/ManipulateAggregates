.. ManipulateAggregates documentation master file, created by
   sphinx-quickstart on Thu Mar 26 12:47:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ManipulateAggregates's documentation!
================================================

.. toctree::
    :maxdepth: 3
    :caption: Contents:

Brief Overview
==============

You are viewing the documentation of ManipulateAggregates, a set of tools for
computational chemistry.
This package comes with two main tools:

1. ``energyscan``:
  - Determines geometricaly diverse local energy minima on the potential energy
    surface of aggregates bound by van-der-Waals or Coulomb forces using
    empirical force fields.

  - Simply run ``energyscan --porphin-example`` or ``energyscan --urea-example`` to
    run the extensive porphin or urea computations published in the following
    paper in your current directory:
    *"A Program for Automatically Predicting Supramolecular Aggregates and
    Its Application to Urea and Porphin"* by Sachse *et al*, accessible at
    https://dx.doi.org/10.1002/jcc.25151

  - You can also run ``energyscan --anthracene-example`` for a quick and less
    demanding scan using the anthracene molecule.

  - Running ``energyscan --longhelp`` will output a complete config file
    including explanations to stdout.

  - If you want to use multiprocessing, set the environment variable
    ``OMP_NUM_THREADS`` to the number of processes you want to use. Happy
    scanning!

2. ``manipagg``:
  - Manipulates internal degrees of freedoms of molecules and aggregates from
    the command line.

  - Computes the electrostatic potential on van-der-Waals surfaces or
    isosurfaces of the electron density based on empirical force fields or
    quantum chemical computations.

  - Simply run ``manipagg --example-vdw`` or ``manipagg --example-iso`` to run an
    example visualization of the electrostatic potential on a molecule's
    van-der-Waals or electrond ensity iso surface, the former as publised in
    the paper *"Introducing double polar heads to highly fluorescent Thiazoles:
    Influence on supramolecular structures and photonic properties"* by
    Kaufmann *et al*, accessible at https://doi.org/10.1016/j.jcis.2018.04.105

  - If you want to use multiprocessing, set the environment variable
    ``OMP_NUM_THREADS`` to the number of processes you want to use. Happy
    rendering and manipulating!


Introduction
============

What you are currently viewing contains the documentation for
ManipulateAggregates, a set of Python scripts written to perform some tasks
related to what I did during my time as a PhD student that will be detailed in
this documentation. For any license-related information, please see the file
called COPYING and the header of each individual ``*.py`` file.

The module ManipulateAggregates consists of four submodules and three
executables (``manipagg``, ``energyscan``, and ``hashsort``). The following is a list
of all all four submodules (in alphabetical order) including a synopsis of the
functionality they provide.

1. collection:
    - read and write several file formats used in computational chemistry:
        - geometry: cube, molden, xyz
        - orbital data: molden
        - volumetric data: cube, dx, xyz
        - frequencies: aims, terachem
        - easily accessible via ``manipagg``

    - control OpenGL from Python:
        - draw (coloured) spheres and trimeshes
        - export to pov-file to render with PoVRay
        - easily accessible via ``manipagg``

    - prefix file names by auto-generated hashes to limit the number of files
      per directory (huge speed-ups for some file systems)

        - used by ``energyscan`` and ``hashsort``

    - read section-less Python config files (with defaults and type checks)

        - used by ``energyscan``

    - most file types can also be read and written when gzipped
    - control gnuplot from Python

2. energyscan:
    - easily accessible via ``energyscan``
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

    - requires `MaAg-bel`_ 
    - requires `FireDeamon`_

3. aggregate:
    - easily accessible via ``manipagg``
    - manipulate (internal) degrees of freedom of molecules and aggregates
    - (compute and) visualize distributions of electrostatic potentials and
      electron densities

        - empirical methods supported via `MaAg-bel`_
        - methods based on results from ab-initio computations supported via
          `FireDeamon`_

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
        - easily accessible via ``manipagg``

    - estimate the character of an orbital by comparison with orbitals of known
      character

    - much of the functionality provided within the submodule orbitalcharacter
      can be sped up when using the C++-library `FireDeamon`_

Example image for the electrostatic potential on a van-der-Waals surface,
quantum chemical data obtained with `TeraChem`_ and rendered by `PoVRay`_:

.. figure:: /images/example_vdw_pot.png


Prerequisites
=============

You need the following dependencies to use this package:
    * Python:
        * `Anaconda`_-based virtual environments highly recommended
        * including the pip package manager
        * including the packages setuptools and wheel
        * tested with version 3.8 but older versions likely work
        * Python 2 working as of 2020, but no efforts are undertaken to stay
          compatible

If you want all the functionality specified above, you also need:
    * `MaAg-bel`_ (modified version of Open Babel)
    * `FireDeamon`_

They in turn have some more dependencies, namely:
    * OpenGL
    * NumPy
    * CGAL
    * GMP
    * MPFR
    * SWIG
    * Eigen3
    * pthreads
    * OpenGL development packages
    * a somewhat recent C++ compiler (C++ 14 must be supported)

There are more completely optional dependencies:
    * `PoVRay`_


Installation
============

With Anaconda on Ubuntu, you can easily install them by doing the following:

.. code-block:: bash

    # Install system packages
    sudo apt-get install libcgal-dev libmpfr-dev libgmp-dev freeglut3 libglu1-mesa-dev
    # If you want to render using PoVRay, run (this is entirely optional):
    sudo apt-get install povray
    # Activate your conda environment
    conda activate <environment>
    # You could also install and activate a new environment like this:
    #   conda create -n manipagg python=3 numpy swig eigen pyopengl \
    #   && conda activate manipagg
    # Install conda packages
    conda install numpy swig eigen pyopengl
    # Install ManipulateAggregates and its dependencies
    pip install ManipulateAggregates

Some environment variables can modify the installation process of `FireDeamon`_.
By default, everything is installed. Important environment variables are:

    * ``MAINST_EIGEN3_DIR`` : set to a non-standard path to the Eigen3 include
      directory. If not set, the setup tries to find the path automatically.

    * ``MAINST_EIGEN3_PREFER_SYSTEM`` : when not specifying ``MAINST_EIGEN3_DIR``,
      set this to ``1`` in order to prefer a system-wide installation of Eigen3
      over a conda installtion, which is preferred otherwise. There is no
      guarantee that a system-wide installation will be found and the conda
      installation should be preferred.

    * ``FDINST_FULL_VIS_SUPPORT`` : if this variable is not ``1`` (``1`` is the
      default), visualization will employ the pyopengl package to interact with
      OpenGL, which is only imported on demand. Effectively, this removes
      ``libGL.so`` as a hard dependency. Thus, you no longer need the system
      packages ``freeglut3`` and ``libglu1-mesa-dev`` and the conda package
      ``pyopengl``.

    * ``FDINST_FULL_SURFACE_SUPPORT`` : if this variable is not ``1`` (``1`` is the
      default), the bare minimum required to run the ``energyscan`` tool will be
      installed. Parts of ``manipagg`` might also work but there are no
      guarantees since ``manipagg`` has a lot more dependencies than ``energyscan``.
      This effectively removes CGAL as a dependency and thus the system packages
      ``libcgal-dev``, ``libmpfr-dev``, and ``libgmp-dev`` are no longer needed as
      well as the conda package ``eigen``.

.. _Anaconda: https://www.anaconda.com/
.. _MaAg-bel: https://github.com/razziel89/MaAg-bel
.. _FireDeamon: https://github.com/razziel89/libfiredeamon
.. _PovRay: http://www.povray.org/
.. _TeraChem: http://petachem.com

If you have all the aforementioned dependencies installed and you environment
activated, simply run:

.. code-block:: bash

    pip install ManipulateAggregates


The `manipagg` Program
======================

.. include:: manipagg.rst


The `energyscan` Program
========================

.. include:: energyscan.rst


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
- License: GNU General Public License version 3
