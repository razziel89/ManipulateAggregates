"""
This module provides wrapper functions to common system calls that allow
the use of files that are in subdirectories according to hashes of their
names.
"""

_open = open
import os, re
from os.path import exists as _exists
import hashlib
from shutil import move as _move

global _DEPTH, _WIDTH, _WIDTHRE, _HASHALG, _HASHFUNC, _HASHES, _LENGTHS

_HASHES = {
          "md5"    : hashlib.md5,
          "sha1"   : hashlib.sha1,
          "sha224" : hashlib.sha224,
          "sha256" : hashlib.sha256,
          "sha384" : hashlib.sha384,
          "sha512" : hashlib.sha512
         }
_LENGTHS = {
          "md5"    : 32,
          "sha1"   : 40,
          "sha224" : 56,
          "sha256" : 64,
          "sha384" : 96,
          "sha512" : 128
         }

_DEPTH=2
_WIDTH=2
_WIDTHRE=re.compile(".{%d}"%(_WIDTH))
_HASHALG="md5"
_HASHFUNC=hashlib.md5

def set_depth(d):
    global _DEPTH
    if d<0:
        raise ValueError("DEPTH must be >0")
    elif d==0:
        _DEPTH=0
        return
    if _LENGTHS[_HASHALG]/d<_DEPTH:
        raise ValueError("Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."%(
            _HASHALG,d,_WIDTH))
    _DEPTH=d
def set_width(w):
    global _WIDTHRE,_WIDTH
    if w<0:
        raise ValueError("WIDTH must be >0")
    elif w==0:
        _WIDTH=0
        return
    if _LENGTHS[_HASHALG]/w<_DEPTH:
        raise ValueError("Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."%(
            _HASHALG,_DEPTH,w))
    _WIDTHRE=re.compile(".{%d}"%(w))
    _WIDTH=w
def set_hashalg(name):
    global _HASHALG,_HASHFUNC
    if _DEPTH==0 or _WIDTH==0:
        return
    if _LENGTHS[name]/_WIDTH<_DEPTH:
        raise ValueError("Chosen hash-function %s does not create hashes long enough for a depth of %d with a directory name width of %d."%(
            name,_DEPTH,_WIDTH))
    _HASHFUNC=_HASHES[name]
    _HASHALG=name

def exists(pathname,nulldepth=False):
    """
    Allow for checking whether a file exists whose name shall be hashed.

    pathname: string, the path to the file to be checked
    nulldepth: bool (optional, default: False), whether or not a depth of
               0 shall be checked and if the file is found there, it shall be moved
               to the appropriate place as if it had been storen in hashed directories.
    """
    if _WIDTH==0 or _DEPTH==0:
        nulldepth=True
    filename    = pathname.split(os.sep)[-1]
    dirname     = pathname[0:-len(filename)-1]
    if not dirname.endswith(os.sep) and len(dirname)>0:
        dirname += os.sep
    hashed      = _HASHFUNC(filename).hexdigest()
    hasheddir   = dirname+os.sep.join([m.string[m.start():m.end()] for m in re.finditer(_WIDTHRE,hashed)][0:_DEPTH])
    hashedcheck = hasheddir+os.sep+filename
    if nulldepth:
        check = dirname+os.sep+filename
        if _exists(check):
            if _exists(hashedcheck):
                raise IOError("Cannot move file, destination already exists: %s -> %s"%(check,hashedcheck))
            if not _exists(hasheddir):
                os.makedirs(hasheddir)
            elif not os.path.isdir(hasheddir):
                raise IOError("File %d already exists but is no directory."%(hasheddir))
            _move(check,hashedcheck)
    return _exists(hashedcheck)

def _listfiles_recursive(dir,depth):
    dirs = [dir]
    for _depth in xrange(depth,0,-1):
        dirs = [d+os.sep+f for d in dirs for f in os.listdir(d) if os.path.isdir(d+os.sep+f)]
        if len(dirs)==0:
            raise IOError("Not enough directories for a depth of %d."%(depth-_depth))
    return [d+os.sep+f for d in dirs for f in os.listdir(d) if os.path.isfile(d+os.sep+f)]

def listfiles(dirname,regex=None,nullsize=True,nulldepth=False):
    """
    List files in a directory and, optinally, only return those files mathing a
    regex, that are not empty or use hashed directories.

    dirname: str, name of the directory
    regex: regular expression (string or regex object), optional (default: None)
           If not None, return only those files whose paths match the given regex.
    nullsize: bool, optional (default: True)
              Whether or not to return objects of zero size.
    nulldepth: bool, optional (default: False)
               Whether or not to ignore hashing.
    """
    if _WIDTH==0 or _DEPTH==0:
        nulldepth=True
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
        return [dirname+os.sep+f for f in os.listdir(dirname)
                if os.path.isfile(dirname+os.sep+f) and regexfunc(f.split(os.sep)[-1]) and nullsizefunc(dirname+os.sep+f)]
    else:
        hashfunc = lambda o: (_HASHFUNC(o.split(os.sep)[-1]).hexdigest()[0:_WIDTH*_DEPTH] == "".join(o.split(os.sep)[-1-_DEPTH:-1]))
        return [f for f in _listfiles_recursive(dirname,_DEPTH) if regexfunc(f.split(os.sep)[-1]) and nullsizefunc(f) and hashfunc(f)]

def hashpath(pathname):
    """
    Return the name of the file prepended by the hashed directory names.

    filename: string, the files name
    """
    if _WIDTH==0 or _DEPTH==0:
        return pathname
    filename    = pathname.split(os.sep)[-1]
    dirname     = pathname[0:-len(filename)-1]
    if not dirname.endswith(os.sep) and len(dirname)>0:
        dirname += os.sep
    hashed      = _HASHFUNC(filename).hexdigest()
    hasheddir   = os.sep.join([m.string[m.start():m.end()] for m in re.finditer(_WIDTHRE,hashed)][0:_DEPTH])
    hashedcheck = dirname+hasheddir+os.sep+filename
    return hashedcheck
