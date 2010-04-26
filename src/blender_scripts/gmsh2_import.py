#!BPY
# -*- coding: utf-8 -*-
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
# Gmsh2 import script for Blender.

"""
Name: 'Gmsh2 (*.msh)'
Blender: 249
Group: 'Import'
Tooltip: 'Import to Gmsh 2 ASCII format'
"""

import Blender
import bpy
import BPyAddMesh

def readGmsh2(filename):
  Blender.Window.WaitCursor(1)
  #mesh = Blender.NMesh.New()
  Vector = Blender.Mathutils.Vector
  
  file = open(filename, 'r')
  
  line = file.readline()
  line = file.readline()
  line = file.readline()
  line = file.readline()
  line = file.readline()
  words = line.split()
  if len(words) != 1:
    Blender.Draw.PupMenu('Error%t|File format error 1')
    return
  num_nodes = int(float(words[0]))
  verts = []
  for i in range(0, num_nodes):
    line = file.readline()
    words = line.split()
    if len(words) != 4:
      Blender.Draw.PupMenu('Error%t|File format error 2')
      return
    x = float(words[1])
    y = float(words[2])
    z = float(words[3])
    verts.append( Vector(x,y,z) )
    
  line = file.readline()
  line = file.readline()
  line = file.readline()
  words = line.split()
  if len(words) != 1:
    Blender.Draw.PupMenu('Error%t|File format error 3')
    return
  num_elements = int(float(words[0]))
  faces = []
  for i in range(0, num_elements):
    line = file.readline()
    words = line.split()
    if len(words) < 3:
      Blender.Draw.PupMenu('Error%t|File format error 4')
      return
    element_type = int(float(words[1]))
    num_tags = int(float(words[2]))
    read_element = 0
    if element_type == 2:
      if len(words) != 3 + 3 + num_tags:
        Blender.Draw.PupMenu('Error%t|File format error 5')
        return
      read_element = 1
    if element_type == 3:
      if len(words) != 3 + 4 + num_tags:
        Blender.Draw.PupMenu('Error%t|File format error 6')
        return
      read_element = 1
    if read_element == 1:
      face_verts = []
      for j in range(3+num_tags, len(words)):
        idx = int(float(words[j])) - 1
        face_verts.append(idx)
      faces.append(face_verts)
      
  BPyAddMesh.add_mesh_simple('gmsh_object', verts, [], faces)
	        
  Blender.Window.WaitCursor(0)
  Blender.Window.RedrawAll()

Blender.Window.FileSelector(readGmsh2, "Import", Blender.sys.makename(ext='.msh'))
