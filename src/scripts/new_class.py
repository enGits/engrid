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
import sys

def get_path(command):
  parts = command.split("/")
  path = ""
  for i in range(len(parts) - 1):
    path = path + parts[i] + "/"
  return path

def get_lic():
  file_name = get_path(sys.argv[0]) + "licence_header.txt"
  read_buf = open(file_name, "r").read()
  lines = read_buf.split("\n")
  lic = ""
  for line in lines:
    lic = lic + "// " + line + "\n"
  return lic
  
def write_h(class_name):
  f = open(class_name.lower() + ".h", "w")
  f.write(get_lic())
  f.write("#ifndef " + class_name.upper() + "_H\n")
  f.write("#define " + class_name.upper() + "_H\n")
  f.write("\n")
  f.write("#include \"engitscloudlib.h\"\n")
  f.write("\n")
  f.write("namespace ECL_NAMESPACE\n")
  f.write("{\n")
  f.write("class " + class_name + ";\n")
  f.write("}\n")
  f.write("\n")
  f.write("namespace ECL_NAMESPACE\n")
  f.write("{\n")
  f.write("\n")
  f.write("class " + class_name + "\n")
  f.write("{\n")
  f.write("};\n")
  f.write("\n")
  f.write("}\n // namespace")
  f.write("\n")
  f.write("#endif // " + class_name.upper() + "_H\n")
  
def write_cpp(class_name):
  f = open(class_name.lower() + ".cpp", "w")
  f.write(get_lic())
  f.write("#include \"" + class_name.lower() + ".h\"\n")
  f.write("\n")
  f.write("namespace ECL_NAMESPACE\n")
  f.write("{\n")
  f.write("\n")
  f.write("} // namespace\n")
  f.write("\n")
  

class_name = sys.argv[1]
write_h(class_name)
write_cpp(class_name)
