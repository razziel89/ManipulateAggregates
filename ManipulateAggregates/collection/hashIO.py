"""Automaticaly generate subdirectories based on file names.

This module provides wrapper functions to common system calls that allow the
use of files that are in subdirectories according to hashes of their names.
This is useful if otherwise many thousands of files would be in the same
directory as some file systems cannot handle such a case.

Imagine a file with the name of "66965_out.dx". Its md5-hash is
1d1696d52e687ade47040221a36673dd. If the hash width has been set to, say, 4
using set_width, sets of 4 letters would be used as subdirectory names.
Using a hash depeth (set using set_depth) of 3, the function hashpath
would return "1d16/96d5/2e68/66965_out.dx" when given the file name.
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
import os, re
from os.path import exists as _exists
import hashlib
from shutil import move as _move

from .p2p3IO import hashstring as _hashstring

global _DEPTH, _WIDTH, _WIDTHRE, _HASHALG, _HASHFUNC, _HASHES, _LENGTHS

_HASHES = {
    "md5": hashlib.md5,
    "sha1": hashlib.sha1,
    "sha224": hashlib.sha224,
    "sha256": hashlib.sha256,
    "sha384": hashlib.sha384,
    "sha512": hashlib.sha512,
}
_LENGTHS = {
    "md5": 32,
    "sha1": 40,
    "sha224": 56,
    "sha256": 64,
    "sha384": 96,
    "sha512": 128,
}

_DEPTH = 2
_WIDTH = 2
_WIDTHRE = re.compile(".{%d}" % (_WIDTH))
_HASHALG = "md5"
_HASHFUNC = hashlib.md5


def set_depth(d):
    """Set the number of levels of hashed directories.

    Args:
        d: (int) number of levels of hashing
    """
    global _DEPTH
    if d < 0:
        raise ValueError("DEPTH must be >0")
    elif d == 0:
        _DEPTH = 0
        return
    if _LENGTHS[_HASHALG] / d < _DEPTH:
        raise ValueError(
            "Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."
            % (_HASHALG, d, _WIDTH)
        )
    _DEPTH = d


def set_width(w):
    """Set the number of levels of hashed directories.

    Args:
        w: (int) number of levels of hashing
    """
    global _WIDTHRE, _WIDTH
    if w < 0:
        raise ValueError("WIDTH must be >0")
    elif w == 0:
        _WIDTH = 0
        return
    if _LENGTHS[_HASHALG] / w < _DEPTH:
        raise ValueError(
            "Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."
            % (_HASHALG, _DEPTH, w)
        )
    _WIDTHRE = re.compile(".{%d}" % (w))
    _WIDTH = w


def set_hashalg(name):
    """Declare the algorithm used to create the hashes.

    You can select from: md5, sha1, sha224, sha256, sha384, sha512
    """
    global _HASHALG, _HASHFUNC
    if _DEPTH == 0 or _WIDTH == 0:
        return
    if _LENGTHS[name] / _WIDTH < _DEPTH:
        raise ValueError(
            "Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."
            % (name, _DEPTH, _WIDTH)
        )
    _HASHFUNC = _HASHES[name]
    _HASHALG = name


def exists(pathname, nulldepth=False):
    """Allow for checking whether a file exists whose name shall be hashed.

    Args:
        pathname: (string) the path to the file to be checked assuming no
            hashing is done.
        nulldepth: (bool) whether or not a depth of 0 shall be checked and if
            the file is found there, it shall be moved to the appropriate place
            as if it had been stored in hashed directories.

    Returns:
        (bool) whether or not the file exists
    """
    if _WIDTH == 0 or _DEPTH == 0:
        nulldepth = True
    filename = pathname.split(os.sep)[-1]
    dirname = pathname[0 : -len(filename) - 1]
    if not dirname.endswith(os.sep) and len(dirname) > 0:
        dirname += os.sep
    hashed = _HASHFUNC(_hashstring(filename)).hexdigest()
    hasheddir = dirname + os.sep.join(
        [m.string[m.start() : m.end()] for m in re.finditer(_WIDTHRE, hashed)][0:_DEPTH]
    )
    hashedcheck = hasheddir + os.sep + filename
    if nulldepth:
        check = dirname + os.sep + filename
        if _exists(check):
            if _exists(hashedcheck):
                raise IOError(
                    "Cannot move file, destination already exists: %s -> %s"
                    % (check, hashedcheck)
                )
            if not _exists(hasheddir):
                os.makedirs(hasheddir)
            elif not os.path.isdir(hasheddir):
                raise IOError(
                    "File %d already exists but is no directory." % (hasheddir)
                )
            _move(check, hashedcheck)
    return _exists(hashedcheck)


def _listfiles_recursive(dir, depth):
    dirs = [dir]
    for _depth in range(depth, 0, -1):
        dirs = [
            d + os.sep + f
            for d in dirs
            for f in os.listdir(d)
            if os.path.isdir(d + os.sep + f)
        ]
        if len(dirs) == 0:
            raise IOError(
                "Not enough directories for a depth of %d." % (depth - _depth)
            )
    return [
        d + os.sep + f
        for d in dirs
        for f in os.listdir(d)
        if os.path.isfile(d + os.sep + f)
    ]


def listfiles(dirname, regex=None, nullsize=True, nulldepth=False):
    """List files in a directory allowing for hashing.

    Optinally, only return those files mathing a regex, those that are not
    empty, or those that use a hash depth of 0.

    Args:
        dirname: (string) name of the directory to be listed
        regex: (string or regex object) if not None, return only those files
            whose paths match the given regex.
        nullsize: (bool) whether or not to return objects of zero size
        nulldepth: (bool) whether or not to ignore hashing (i.e., have depth of 0)

    Returns:
        a list of filenames
    """
    if _WIDTH == 0 or _DEPTH == 0:
        nulldepth = True
    if regex is None:
        regexfunc = lambda o: True
    else:
        regex = re.compile(regex)
        regexfunc = lambda o: (regex.match(o) is not None)
    if nullsize:
        nullsizefunc = lambda o: True
    else:
        nullsizefunc = lambda o: (os.stat(o).st_size > 0)
    if nulldepth:
        return [
            dirname + os.sep + f
            for f in os.listdir(dirname)
            if os.path.isfile(dirname + os.sep + f)
            and regexfunc(f.split(os.sep)[-1])
            and nullsizefunc(dirname + os.sep + f)
        ]
    else:
        hashfunc = lambda o: (
            _HASHFUNC(_hashstring(o.split(os.sep)[-1])).hexdigest()[0 : _WIDTH * _DEPTH]
            == "".join(o.split(os.sep)[-1 - _DEPTH : -1])
        )
        return [
            f
            for f in _listfiles_recursive(dirname, _DEPTH)
            if regexfunc(f.split(os.sep)[-1]) and nullsizefunc(f) and hashfunc(f)
        ]


def hashpath(pathname):
    """Return the name of the file prepended by the hashed directory names.

    Args:
        pathname: (string) the files name whose name shall be hashed

    Returns:
        the hashed path
    """
    if _WIDTH == 0 or _DEPTH == 0:
        return pathname
    filename = pathname.split(os.sep)[-1]
    dirname = pathname[0 : -len(filename) - 1]
    if not dirname.endswith(os.sep) and len(dirname) > 0:
        dirname += os.sep
    hashed = _HASHFUNC(_hashstring(filename)).hexdigest()
    hasheddir = os.sep.join(
        [m.string[m.start() : m.end()] for m in re.finditer(_WIDTHRE, hashed)][0:_DEPTH]
    )
    hashedcheck = dirname + hasheddir + os.sep + filename
    return hashedcheck
