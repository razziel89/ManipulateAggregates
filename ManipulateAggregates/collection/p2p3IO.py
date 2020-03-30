"""A consistent collection of wrapper functions for opening, closing and writing to files using Python 2 and 3.

Files are limited to ASCII encoding as other encodings are generally not supported by ManipulateAggregates.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import io

__ENCODING__ = "ascii"
_open = open
# start feature detection
try:
    # Python2
    unicode("1")
    tounicode = lambda s: unicode(s)
    csvwriteopen = lambda path: _open(path, "wb")
except NameError:
    # Python3
    tounicode = lambda s: str(s)
    csvwriteopen = lambda path: _open(path, "w", encoding=__ENCODING__, newline="")
try:
    # Python3
    bytes("asdf", "ascii")
    hashstring = lambda s: bytes(s, __ENCODING__)
    tobasestring = lambda s: s.decode(encoding=__ENCODING__)
except TypeError:
    # Python2
    hashstring = lambda s: str(s)
    tobasestring = lambda s: s
try:
    # Python2
    isinstance("123", basestring)
    isbasestring = lambda s: isinstance(s, basestring)
except NameError:
    # Python3
    isbasestring = lambda s: isinstance(s, str)


def open(path, mode="r"):
    """Open file located at path in mode mode (always text mode).
    
    Binary mode is not supported.
    """
    if "b" in mode:
        raise ValueError("This 'open' function does not support binary mode.")
    return io.open(path, mode, encoding=__ENCODING__)


def writeto(handle, text):
    """Writes test to handle."""
    handle.write(tounicode(text))


def close(handle):
    """Close handle."""
    handle.close()
