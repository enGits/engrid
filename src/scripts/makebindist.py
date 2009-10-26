#!/usr/bin/env python

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
