"""A consistent collection of wrapper functions for opening, closing and writing to files using Python 2 and 3.

Files are limited to ASCII encoding as other encodings are generally not supported by ManipulateAggregates.

@package ManipulateAggregates.aggregate.p2p3IO
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import io
##\cond
__ENCODING__="ascii"
try:
    unicode("1")
except NameError:
    unicode=lambda s:str(s)
##\endcond

def open(path,mode="r"):
    """Open file located at @a path in @a mode mode (always text mode).
    
    Binary mode is not supported.
    """
    if "b" in mode:
        raise ValueError("This 'open' function does not support binary mode.")
    io.open(path,mode,encoding=__ENCODING__)

def writeto(handle,text):
    """Writes @a test to @a handle."""
    handle.write(unicode(text))

def close(handle):
    """Close @a handle."""
    handle.close()
