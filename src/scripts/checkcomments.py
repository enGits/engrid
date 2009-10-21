#!/usr/bin/env python
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008,2009 Oliver Gloth                                     +
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
# + along with enGrid. If not, see <http:#www.gnu.org/licenses/>.        +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# DESCRIPTION:
# This script outputs all lines prefixed with "///@@@" in the input files given.

import sys

def trimline(line):
  append = False
  trimmedline = ''
  for i in range (0,len(line)):
    #print line[i]
    if line[i] != ' ':
      append = True
    if append:
      trimmedline = trimmedline + line[i]
  return trimmedline
    
def fileinfo(name,ln):
  info = str(ln)
  while len(info) < 5:
    info = info + ' '
  info = name + ' line ' + info;
  return info    

for i in range(1,len(sys.argv)):
  f = open(sys.argv[i])
  #print f
  line = f.readline()
  conti = True
  ln = 1
  while line:
    tline = trimline(line)
    #print tline
    if tline[0:6] == '///@@@':
      print fileinfo(sys.argv[i],ln) + ':' + tline[6:len(tline)-1]
      conti = True
    else:
      conti = False
    line = f.readline()
    ln = ln + 1
  f.close
