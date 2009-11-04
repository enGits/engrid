#!BPY

"""
Name: 'Gmsh2 (*.msh)'
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

