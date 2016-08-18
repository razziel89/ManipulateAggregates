#!/bin/env python
import sys,os,re
from collection import hashIO

if __name__ == "__main__":
    if len(sys.argv)!=3:
        raise ValueError("Not enough arguments given, must be 2.")
    dir = sys.argv[1]
    if not os.path.isdir(dir):
        raise IOError("First argument must be a directory")
    dxregex = re.compile("^[1-9][0-9]*_%s$"%(sys.argv[2]))
    count = 0
    for f in os.listdir(dir):
        if os.path.isfile(dir+os.sep+f) and dxregex.match(f) is not None:
            count += 1
            print "%s -> %s"%(dir+os.sep+f,hashIO.hashpath(dir+os.sep+f))
            hashIO.exists(dir+os.sep+f,nulldepth=True)
    print "Sorted %d files into hashed directories."%(count)
