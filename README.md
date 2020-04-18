# Overview

You are viewing the readme of ManipulateAggregates, a set of tools for
computational chemistry.
This package comes with two main tools:

* `energyscan`:
  * Determines geometricaly diverse local energy minima on the potential energy
    surface of aggregates bound by van-der-Waals or Coulomb forces using
    empirical force fields.
  * Simply run `energyscan --porphin-example` or `energyscan --urea-example` to
    run the extensive porphin or urea computations published in the following
    paper in your current directory:
    *"A Program for Automatically Predicting Supramolecular Aggregates and
    Its Application to Urea and Porphin"* by Sachse *et al*, accessible at
    https://dx.doi.org/10.1002/jcc.25151
  * You can also run `energyscan --anthracene-example` for a quick and less
    demanding scan using the anthracene molecule.
  * Running `energyscan --longhelp` will output a complete config file
    including explanations to stdout.
  * If you want to use multiprocessing, set the environment variable
    `OMP_NUM_THREADS` to the number of processes you want to use. Happy
    scanning!

* `manipagg`:
  * Manipulates internal degrees of freedoms of molecules and aggregates from
    the command line.
  * Computes the electrostatic potential on van-der-Waals surfaces or
    isosurfaces of the electron density based on empirical force fields or
    quantum chemical computations.
  * Simply run `manipagg --example-vdw` or `manipagg --example-iso` to run an
    example visualization of the electrostatic potential on a molecule's
    van-der-Waals or electrond ensity iso surface, the former as publised in
    the paper *"Introducing double polar heads to highly fluorescent Thiazoles:
    Influence on supramolecular structures and photonic properties"* by
    Kaufmann *et al*, accessible at https://doi.org/10.1016/j.jcis.2018.04.105
  * If you want to use multiprocessing, set the environment variable
    `OMP_NUM_THREADS` to the number of processes you want to use. Happy
    rendering and manipulating!

Please see the documentation for a detailed description and a full list of
features on <https://razziel89.github.io/ManipulateAggregates/> (provided via
GitHub pages).

# Quick installation

If you are running Ubuntu and use Anaconda to manage your Python environments,
you can easily install ManipulateAggregates the following way:

```bash
# Install system packages
sudo apt-get install libcgal-dev libmpfr-dev libgmp-dev freeglut3 libglu1-mesa-dev
# If you want to render using PoVRay, run:
sudo apt-get install povray
# Install and activate a new environment like this:
conda create -n manipagg python=3 numpy swig eigen pyopengl
conda activate manipagg
# Install ManipulateAggregates and its dependencies
pip install ManipulateAggregates
```

Please refer to the documentation on
<https://razziel89.github.io/ManipulateAggregates/#prerequisites> for more
information.

# Contributing

Contributions are very welcome!
Please simply open a pull request.
If you would like to make large-ish contributions, it might be prudent to first
contact the maintainer to better co-ordinate those efforts.

This project uses the following auto-formatter:
* Python code: black (the uncompromising Python code formatter)
  <https://github.com/psf/black>

Please make sure to auto-format your pull request with those options.
Furthermore, please document any code you add.
Happy contributing!
