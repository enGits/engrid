#!/usr/bin/env python
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
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
import os

def update_file(file_name, header, comment, num_preserve):
  #print file_name
  f = open(file_name, "r")
  old_lines = f.readlines()
  f.close()
  new_lines = []
  in_header = True
  for i in range(num_preserve):
    new_lines.append(old_lines[i])
  for i in range(num_preserve, len(old_lines)):
    line = old_lines[i]
    if in_header:
      if not line.startswith(comment):
        in_header = False
    if not in_header:
      new_lines.append(line)
  f = open(file_name, "w")
  for i in range(num_preserve):
    f.write(new_lines[i])
  for line in header:
    f.write(comment + " " + line)
  for i in range(num_preserve, len(new_lines)):
    f.write(new_lines[i])
  

def recursive_traversal(dir, header):
  fns = os.listdir(dir)
  for fn in fns:
    if fn != "tetgen":
      file_name = os.path.join(dir,fn)
      if os.path.isdir(file_name):
        recursive_traversal(file_name, header)
      else:
        if (file_name.endswith(".h")):
          update_file(file_name, header, "//", 0)
        if (file_name.endswith(".cpp")):
          update_file(file_name, header, "//", 0)
        if (file_name.endswith(".cxx")):
          update_file(file_name, header, "//", 0)
        if (file_name.endswith(".cu")):
          update_file(file_name, header, "//", 0)
        if (file_name.endswith(".pro")):
          update_file(file_name, header, "#", 0)
        if (file_name.endswith(".pri")):
          update_file(file_name, header, "#", 0)
        if (file_name.endswith(".sh")):
          update_file(file_name, header, "#", 1)
        if (file_name.endswith(".bash")):
          update_file(file_name, header, "#", 1)
        if (file_name.endswith(".py")):
          update_file(file_name, header, "#", 1)
        if (file_name.endswith(".f")):
          update_file(file_name, header, "c    ", 0)
        if (file_name.endswith("CMakeLists.txt")):
          update_file(file_name, header, "#", 0)


f = open("licence_header.txt")
header = f.readlines()
f.close

recursive_traversal(".", header)


