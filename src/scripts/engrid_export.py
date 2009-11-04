#!BPY

"""
Name: 'Engrid (*.begc)'
Blender: 249
Group: 'Export'
Tooltip: 'Export to Engrid'
"""

import Blender
import bpy

def write(filename):
  if not filename.lower().endswith('.begc'):
    filename += '.begc'
  out = file(filename, "w")
  objects = Blender.Object.GetSelected()
  
  num_objects = 0
  for object in objects:
    if object.type == 'Mesh':
      num_objects = num_objects + 1
      
  out.write('%d\n' % num_objects)
  node_offset = 0
  for object in objects:
    if object.type == 'Mesh':
      out.write(object.name)
      out.write('\n')
  for object in objects:
    if object.type == 'Mesh':
      mesh  = object.getData(0,1)
      faces = mesh.faces
      nodes = mesh.verts
      out.write('%d' % len(nodes))
      out.write(' %d\n' % len(faces))
      for n in nodes:
        out.write("%e " % n.co[0])
        out.write("%e " % n.co[1])
        out.write("%e\n" % n.co[2])
      for f in faces:
        N = len(f.verts)
        if N < 3 and N > 4:
          Blender.Draw.PupMenu('Error%t|Only triangles and quads allowed')
          return
        out.write("%d" % N)
        for v in f.verts:
          out.write(' %d' % (v.index + node_offset))
        out.write('\n')
      node_offset = node_offset + len(nodes)

Blender.Window.FileSelector(write, "Export", Blender.sys.makename(ext='.begc'))

