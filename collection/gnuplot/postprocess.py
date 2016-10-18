"""
This module defines a set of functions that can be applied to data.
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

import numpy as np

def _normalize(data,xcol,ycols,retvals,start=None,stop=None):
    if start is None:
        start = np.min(data[xcol])
    if stop is None:
        stop = np.max(data[xcol])
    start,stop = min((start,stop)),max((start,stop))
    for y in ycols:
        maximum = np.max(data[y])
        if maximum != 0.0:
            data[y] /= maximum
    return data
    
def _normalize_value(data,xcol,ycols,retvals,value):
    left   = sorted(data.T[data[xcol]<=value],key=lambda e:e[xcol])[-1]
    right  = sorted(data.T[data[xcol]>value],key=lambda e:e[xcol])[0]
    leftx  = left[0]
    left   = left[1:]
    rightx = right[0]
    right  = right[1:]
    weight = (value-leftx)/(leftx-rightx)
    scale  = (1.0-weight)*left + weight*right
    for y,s in zip(ycols,scale):
        data[y] /= s
    return data

def _xshift(data,xcol,ycols,retvals,shift):
    pass

def _xscale(data,xcol,ycols,retvals,scale):
    pass

def _yshift(data,xcol,ycols,retvals,shift):
    pass

def _yscale(data,xcol,ycols,retvals,scale):
    pass

def _nothing(data,xcol,retvals,ycols):
    return data

functions = {
        "normalize"         :   _normalize,
        "normalize_value"   :   _normalize_value,
        "xshift"            :   _xshift,
        "yshift"            :   _yshift,
        "xscale"            :   _xscale,
        "yscale"            :   _yscale,
        "nothing"           :   _nothing
        }

def apply(data,postprocess,args,xcol,ycols,ignore_returns=True):
    """
    Apply some postprocessing functions in order to some data and return the
    data and the return values of the functions if so desired. A plain python
    list will always be returned if postprocess has been set but numpy will be
    used in the meantime.

    data: list of lists of floats
        The data that is to be manipulated.
    postprocess: list of strings
        The names of the manipulation routines that shall be applied to the data.
    args: list of lists of appropriate data
        The arguments that shall be passed to the manipulation routines.
    xcol: int
        The coloumn that shall be considered the argument to which function values
        are asigned.
    ycols: list of ints
        The coloumns that contain the function values.
    ignore_returns: bool, optional, default: True
        Whether or not to ignore return values from the manipulation routines.
        This is usually only needed for debugging purposes.
    """
    retvals = []
    np_data = np.array(data).T
    try:
        for p,a in zip(postprocess,args):
            np_data = functions[p](np_data,xcol,ycols,retvals,*a)
    except KeyError as e:
        raise KeyError("Unknown function \"%s\" for postprocessing data."%(p),e)
    if not ignore_returns:
        return np_data.T, retvals
    else:
        return np_data.T
