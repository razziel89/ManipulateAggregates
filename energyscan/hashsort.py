#!/bin/env python
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
