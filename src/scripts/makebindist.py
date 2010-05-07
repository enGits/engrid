#!/usr/bin/env python
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2010 enGits GmbH                                     +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 

import os
import tempfile
import re
import shutil

#filename = tempfile.mkstemp()[1]
tmpdir = tempfile.mkdtemp()

print "================"
print tmpdir
print "================"

ldd_log = tmpdir + "/ldd.out"
command = "ldd ./engrid > " + ldd_log
print command
os.system(command)

libdir=tmpdir + "/enGrid"
os.mkdir(libdir)

count = 0

liblist = []

infile = open(ldd_log, "r")
if infile:
    print "name="+infile.name
    for line in infile:
        p1 = re.compile('\s\S+\s=>\s(\S+)\s\(\S+\)\s')
        p2 = re.compile('\s(\S+)\s\(\S+\)\s')

        m1 = p1.match(line)
        if m1:
            print m1.group(1)
            count += 1
            liblist.append(m1.group(1))
            
        m2 = p2.match(line)
        if m2:
            print m2.group(1)
            count += 1
            liblist.append(m2.group(1))

    infile.close()

print "count = ", count
print "liblist = ", liblist
print "len(liblist) = ", len(liblist)

for lib in liblist:
  print lib
  shutil.copy(lib,libdir)

print "================"
print tmpdir
print "================"

#cp engrid $TMPDIR/enGrid
#cp start_engrid enGrid
#tar -czvf enGrid_bin.tar.gz $TMPDIR/enGrid/*
# rm -rf enGrid
# rm ldd.out
