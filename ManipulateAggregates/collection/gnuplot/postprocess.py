"""This module defines a set of functions that can be applied to data.
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
import logging

logger = logging.getLogger(__name__)
try:
    import numpy
except ImportError:
    logger.warning("Could not import numpy")

global functions


def apply(data, postprocess, args, xcol, ycols, ignore_returns=True):
    """Main function to postprocess data that was converted to gnuplot format.

    Apply some postprocessing functions in order to some data and return the
    data and the return values of the functions if so desired. A plain python
    list will always be returned if postprocess has been set but numpy will be
    used in the meantime.

    Args:
        data: (list of lists of floats) the data that is to be manipulated.
        postprocess: (list of strings) the names of the manipulation routines that
            shall be applied to data
        args: (list of lists of appropriate data) the arguments that shall be
            passed to the manipulation routines. See their own documentations.
        xcol: (int) the coloumn that shall be considered the argument to which
            function values are asigned
        ycols: (list of ints) the coloumns that contain the function values
        ignore_returns: (bool) whether or not to ignore return values from the
            manipulation routines. This is usually only needed for debugging purposes.
    """
    global functions
    retvals = []
    np_data = numpy.array(data).T
    try:
        for p, a in zip(postprocess, args):
            np_data = functions[p](np_data, xcol, ycols, retvals, *a)
    except KeyError as e:
        raise KeyError('Unknown function "%s" for postprocessing data.' % (p), e)
    if not ignore_returns:
        return np_data.T, retvals
    else:
        return np_data.T


def normalize(data, xcol, ycols, retvals, start=None, stop=None):
    """Normalize the data.

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        start: (float) the start of the normalization range
        stop: (float) the end of the normalization range

    """
    if start is None:
        start = numpy.min(data[xcol])
    if stop is None:
        stop = numpy.max(data[xcol])
    start, stop = min((start, stop)), max((start, stop))
    for y in ycols:
        maximum = numpy.max(data[y])
        if maximum != 0.0:
            data[y] /= maximum
    return data


def normalize_value(data, xcol, ycols, retvals, value):
    """Scale the data so that the value associated with a certain x position is 1.

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        value: (float) the position that shall have the value 1

    """
    left = sorted(data.T[data[xcol] <= value], key=lambda e: e[xcol])[-1]
    right = sorted(data.T[data[xcol] > value], key=lambda e: e[xcol])[0]
    leftx = left[0]
    left = left[1:]
    rightx = right[0]
    right = right[1:]
    weight = (value - leftx) / (leftx - rightx)
    scale = (1.0 - weight) * left + weight * right
    for y, s in zip(ycols, scale):
        data[y] /= s
    return data


def xshift(data, xcol, ycols, retvals, shift):
    """NOT YET IMPLEMENTED

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        shift: NOT YET USED
    """
    pass


def xscale(data, xcol, ycols, retvals, scale):
    """NOT YET IMPLEMENTED

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        scale: NOT YET USED
    """
    pass


def yshift(data, xcol, ycols, retvals, shift):
    """NOT YET IMPLEMENTED

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        shift: NOT YET USED
    """
    pass


def yscale(data, xcol, ycols, retvals, scale):
    """NOT YET IMPLEMENTED

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
        scale: NOT YET USED
    """
    pass


def nothing(data, xcol, ycols, retvals):
    """Do nothing to the data.

    Args:
        data: automatically set by apply
        xcol: automatically set by apply
        ycols: automatically set by apply
        retvals: automatically set by apply
    """
    return data


## a dictionary linking function pointers to postprocessing names
functions = {
    "normalize": normalize,
    "normalize_value": normalize_value,
    "xshift": xshift,
    "yshift": yshift,
    "xscale": xscale,
    "yscale": yscale,
    "nothing": nothing,
}
