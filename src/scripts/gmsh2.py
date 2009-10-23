#!BPY
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
# Export script for Blender.

"""
Name: 'Gmsh2'
Blender: 249
Group: 'Export'
Tooltip: 'Export to Gmsh 2 ASCII format'
"""

import Blender
import bpy

def writeGmsh2(filename):
  if not filename.lower().endswith('.msh'):
    filename += '.msh'
  out = file(filename, "w")
  scn = bpy.data.scenes.active
  object = scn.objects.active
  if not object:
    Blender.Draw.PupMenu('Error%t|Select 1 active object')
    return
  if object.type != 'Mesh':
    Blender.Draw.PupMenu('Error%t|Select a mesh object')
    return
  
  mesh  = object.getData(0,1)
  faces = mesh.faces
  nodes = mesh.verts
  
  out.write('$MeshFormat\n2 0 8\n$EndMeshFormat\n')
  out.write('$Nodes\n%d\n' % len(nodes))
  
  i = 1
  for n in nodes:
    out.write("%d " % i)
    out.write("%e " % n.co[0])
    out.write("%e " % n.co[1])
    out.write("%e\n" % n.co[2])
    i = i + 1
    
  out.write('$EndNodes\n')
  out.write('$Elements\n%d\n' % len(faces))
  i = 1
  for f in faces:
    out.write("%d " % i)
    N = len(f.verts)
    if N < 3 and N > 4:
      Blender.Draw.PupMenu('Error%t|Only triangles and quads allowed')
      return
    if N == 3:
      out.write('2 ')
    if N == 4:
      out.write('3 ')
    out.write('1 %d' % (f.mat + 1))
    for v in f.verts:
      out.write(' %d' % (v.index + 1))
    out.write('\n')    
    i = i + 1
    
  out.write('$EndElements\n')

Blender.Window.FileSelector(writeGmsh2, "Export", Blender.sys.makename(ext='.msh'))

