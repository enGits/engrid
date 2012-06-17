#!BPY
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2011 enGits GmbH                                     +
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

"""
Name: 'Engrid (*.begc)'
Blender: 249
Group: 'Import'
Tooltip: 'Import from Engrid'
"""

import Blender
import bpy
import BPyAddMesh

def readEngrid(filename):
  Blender.Window.WaitCursor(1)
  Vector = Blender.Mathutils.Vector
  
  in_file = file(filename, "r")
  line = in_file.readline()
  Nobjects = int(line)
  #print "Nobjects=",Nobjects
  object_names = []
  for i in range(0, Nobjects):
    line = in_file.readline()
    object_names.append(line.strip())
  
  #print "object_names=",object_names
  
  global_verts = []
  offset = 0
  for i_object in range(0, Nobjects):
    line = in_file.readline()
    words = line.split()
    Nverts = int(words[0])
    Nfaces = int(words[1])
    #print "Nverts=",Nverts
    #print "Nfaces=",Nfaces
    
    local_verts = []
    for i_vert in range(0, Nverts):
      line = in_file.readline()
      words = line.split()
      x = float(words[0])
      y = float(words[1])
      z = float(words[2])
      #print "x,y,z=",x,y,z
      local_verts.append( Vector(x,y,z) )
      global_verts.append( Vector(x,y,z) )
    
    faces = []
    for i_face in range(0, Nfaces):
          line = in_file.readline()
          words = line.split()
          if len(words) < 3:
            Blender.Draw.PupMenu('Error%t|File format error 4')
            return
          Nverts_in_face = int(words[0])
          if len(words) != 1 + Nverts_in_face:
            Blender.Draw.PupMenu('Error%t|File format error 5')
            return
          face_verts = []
          for i_face_vert in range(0, Nverts_in_face):
            idx = int(words[i_face_vert + 1]) - offset
            face_verts.append(idx)
          #print "face_verts=",face_verts
          faces.append(face_verts)
    
    #print "Adding object ",object_names[i_object]
    BPyAddMesh.add_mesh_simple(object_names[i_object], local_verts, [], faces)

    offset += Nverts

  in_file.close()

  Blender.Window.WaitCursor(0)
  Blender.Window.RedrawAll()

Blender.Window.FileSelector(readEngrid, "Import", Blender.sys.makename(ext='.begc'))
