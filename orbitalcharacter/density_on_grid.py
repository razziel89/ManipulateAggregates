"""
Compute the electron density on a grid. First call prepare_grid_calculation and
then call density_on_grid(coefficients_list,data) where data is what prepare_grid_calculation
returned.
"""

import sys
import os
from multiprocessing import Pool, Event

import numpy as np

global data

def prepare_grid_calculation(grid,basis):
    """
    Create data structures suitable for efficiently computing
    the elctron density on an arbitrary grid. Call this first
    and then density_on_grid(coefficients,data) where data
    is what this function returns.

    grid: list of [float,float,float]
        The Cartesian coordinates of the grid
    basis: a list of [A,L,Prim]
           with
           A: a list of 3 floats
                The center of the contracted Cartesian Gaussian function
           L: a list of 3 ints
                The polynomial exponents of the contracted Cartesian Gaussian
           Prim: a list of [alpha,pre]
                with
                alpha: float
                    The exponential factor of the primitive Gaussian function
                pre: float
                    The contraction coefficient of the primitive Gaussian function
    """
    newgrid       = np.array(grid)
    primcoords    = np.array([ b[0] for b in basis for prim in xrange(len(b[2])) ],dtype=float)
    primaxyz      = np.array([ b[1] for b in basis for prim in xrange(len(b[2])) ],dtype=int)
    primalpha     = np.array([ prim[0] for b in basis for prim in b[2] ],dtype=float)
    primcoeffs    = np.array([ prim[1] for b in basis for prim in b[2] ],dtype=float)
    indices       = np.array([ bc for b,bc in zip(basis,xrange(len(basis))) for prim in xrange(len(b[2])) ],dtype=int)
    return newgrid,primcoords,primaxyz,primcoeffs,primalpha,indices

def _density_on_grid_parallel_init(grid,primcoords,primaxyz,primcoeffs,primalpha,indices,terminating):
    global data
    data = (grid,primcoords,primaxyz,primcoeffs,primalpha,indices,terminating)

def _density_on_grid_process(coefficients):
    """
    This function controls all the worker processes that calculate the
    electron density on an arbitrary grid. The formula is quite simple:
    rho(r) = sum of |phi(r)|^2 over all basis functions phi
    """
    global data
    grid,primcoords,primaxyz,primcoeffs,primalpha,indices,terminating = data
    try:
        if not terminating.is_set():
            length = len(primcoeffs)
            gridlength = len(grid)
            print gridlength
            result = [0.0]*gridlength
            count=0
            np_coeff         = np.array(coefficients)[indices]
            negprimexp       = -primalpha
            coorddiff        = np.zeros((length,3),dtype=float)
            coorddiff_2      = np.zeros((length,3),dtype=float)
            absvals_2        = np.zeros((length,),dtype=float)
            absvals_2_mult   = np.zeros((length,),dtype=float)
            exponents        = np.zeros((length,),dtype=float)
            prefactor_single = np.zeros((length,3),dtype=float)
            prefactor        = np.zeros((length,),dtype=float)
            product_1        = np.zeros((length,),dtype=float)
            product_2        = np.zeros((length,),dtype=float)
            for r in grid:
                #START: calculate exp(-alpha*|r-A|^2)
                np.subtract(r,primcoords,out=coorddiff)
                np.square(coorddiff,out=coorddiff_2)
                np.sum(coorddiff_2,out=absvals_2,axis=1)
                np.multiply(negprimexp,absvals_2,out=absvals_2_mult)
                np.exp(absvals_2_mult,out=exponents)
                #END:   calculate exp(-alpha*|r-A|^2)
                #START: calculate (x-X)^ax * (y-Y)^ay * (z-Z)^az
                np.power(coorddiff,primaxyz,out=prefactor_single)
                np.prod(prefactor_single,out=prefactor,axis=1)
                #END:   calculate (x-X)^ax * (y-Y)^ay * (z-Z)^az
                #START: multiply MO-coefficient and contraction-coefficient and prefactor and exponent
                np.multiply(prefactor,exponents,out=product_1)
                np.multiply(product_1,primcoeffs,out=product_2)
                np.multiply(product_2,np_coeff,out=product_1)
                #END:   multiply MO-coefficient and contraction-coefficient and prefactor and exponent
                result[count] = pow(np.sum(product_1),2)
                count+=1
        return result
    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

def density_on_grid(coefficients_list,data):
    """
    Calculate the electron density on an arbitrary grid. Multiprocessing is supported
    for several densities at the same time. The environment variable OMP_NUM_THREADS
    declares how many processes can be spawned at the same time.

    coefficients_list: list of lists of floats
        The expansion coefficients for all the wafefunctions whose density
        is to be computed.
    data: what prepare_grid_calculation returned
    """
    grid,primcoords,primaxyz,primcoeffs,primalpha,indices = data
    try:
        nr_threads = int(os.environ["OMP_NUM_THREADS"])
    except KeyError:
        nr_threads = 1
    if nr_threads > len(coefficients_list):
        nr_threads = len(coefficients_list)
    terminating = Event()
    pool = Pool(nr_threads, initializer=_density_on_grid_parallel_init, initargs=(
        grid,primcoords,primaxyz,primcoeffs,primalpha,indices,terminating))
    result = pool.map(_density_on_grid_process, coefficients_list)
    return result
