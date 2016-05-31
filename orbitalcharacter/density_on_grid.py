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

def prepare_grid_calculation(grid,basis,scale=1.0,type='c++'):
    """
    Create data structures suitable for efficiently computing
    the elctron density on an arbitrary grid. Call this first
    and then density_on_grid(coefficients_list,data) where data
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
    scale: float
        Every coordinate of every gridpoint will be divided by this value.
    type: str
        Either c++ (for the fast c++ variant of the algorithm) or
        numpy (for the older numpy version).
    computation: str
        Either 'density' or 'potential' depending on what shall be computed
    """
    if type=='c++':
        from FireDeamon import InitializeElectronDensityPy
        return InitializeElectronDensityPy(grid,basis,scale=scale)
    elif type=='numpy':
        newgrid       = np.array(grid)/scale
        primcoords    = np.array([ b[0] for b in basis for prim in xrange(len(b[2])) ],dtype=float)
        primaxyz      = np.array([ b[1] for b in basis for prim in xrange(len(b[2])) ],dtype=int)
        primalpha     = np.array([ prim[0] for b in basis for prim in b[2] ],dtype=float)
        primcoeffs    = np.array([ prim[1] for b in basis for prim in b[2] ],dtype=float)
        indices       = np.array([ bc for b,bc in zip(basis,xrange(len(basis))) for prim in xrange(len(b[2])) ],dtype=int)
        return newgrid,primcoords,primaxyz,primcoeffs,primalpha,indices
    else:
        raise ValueError("Wrong type for calculation given.")

def _density_on_grid_parallel_init(grid,primcoords,primaxyz,primcoeffs,primalpha,indices,volume,terminating):
    global data
    data = (grid,primcoords,primaxyz,primcoeffs,primalpha,indices,volume,terminating)

def _density_on_grid_process(coefficients):
    """
    This function controls all the worker processes that calculate the
    electron density on an arbitrary grid. The formula is quite simple:
    rho(r) = sum of |phi(r)|^2 over all basis functions phi
    """
    global data
    grid,primcoords,primaxyz,primcoeffs,primalpha,indices,volume,terminating = data
    try:
        if not terminating.is_set():
            length = len(primcoeffs)
            gridlength = len(grid)
            result = np.zeros((gridlength,),dtype=float)
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
                result[count] = pow(np.sum(product_1),2)/volume
                count+=1
        return result
    except KeyboardInterrupt:
        print >>sys.stderr, "Terminating worker process "+str(os.getpid())+" prematurely."

def density_on_grid(coefficients_list,data,volume=1.0,async=True,type='c++',normalize_to=None, cutoff=-1.0):
    """
    Calculate the electron density on an arbitrary grid. Multiprocessing is supported
    for several densities at the same time. The environment variable OMP_NUM_THREADS
    declares how many processes can be spawned at the same time.

    coefficients_list: list of lists of floats
        The expansion coefficients for all the wafefunctions whose density
        is to be computed.
    data: what prepare_grid_calculation returned
    volume: float
        Scale volumetric data by the inverse of this volume
        to get a true density in the voxel.
    async: boolean
        Return data in arbitrary order (with respect to the entries in
        coefficients_list) but print progress reports. If type is c++, the
        order is preserved even if True.
    type: str
        Either c++ (for the fast c++ variant of the algorithm) or
        numpy (for the older numpy version).
    normalize_to: float
        If not None, make it so that the sum over all returned values
        is equal to the given number.
    cutoff: float (in units of the grid!!!)
        If a point on the grid and the center of a basis function are farther
        apart from each other than this value, the density will not be evaluated.
        Switch-off use of cutoff by setting a negative value (default). Only works
        with the C++ algorithm since it would not give any speedup for the numpy one.
    """
    if type=='c++':
        from FireDeamon import ElectronDensityPy
        result = np.array(ElectronDensityPy(coefficients_list,data,volume=volume,prog_report=async,detailed_prog=False,cutoff=cutoff))
    elif type=='numpy':
        grid,primcoords,primaxyz,primcoeffs,primalpha,indices = data
        try:
            nr_threads = int(os.environ["OMP_NUM_THREADS"])
        except KeyError:
            nr_threads = 1
        if nr_threads > len(coefficients_list):
            nr_threads = len(coefficients_list)
        terminating = Event()
        pool = Pool(nr_threads, initializer=_density_on_grid_parallel_init, initargs=(
            grid,primcoords,primaxyz,primcoeffs,primalpha,indices,volume,terminating))
        try:
            if async:
                result=[]
                rcount=0
                maxrcount=len(coefficients_list)
                print "0/"+str(maxrcount)
                for temp in pool.imap_unordered(_density_on_grid_process, coefficients_list):
                    rcount+=1
                    reportstring = str(rcount)+"/"+str(maxrcount)
                    print reportstring
                    result.append(temp)
            else:
                result = pool.map(_density_on_grid_process, coefficients_list)
            pool.close()
        except KeyboardInterrupt:
            print >>sys.stderr,"Caught keyboard interrupt."
            pool.terminate()
            print >>sys.stderr,"Terminating main routine prematurely."
        finally:
            pool.join()
        result = np.array(result)
        result = np.sum(result,axis=1)
        #result = [sum(r) for r in zip(*result)]
    else:
        raise ValueError("Wrong type for calculation given.")
    if normalize_to is not None:
        result *= normalize_to / np.sum(result)
    return result
