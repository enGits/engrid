# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2013 enGits GmbH                                      +
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

bl_info = {
    "name": "Export to enGrid (.begc)",
    "author": "Oliver Gloth (enGits GmbH)",
    "version": (1, 0),
    "blender": (2, 6, 3),
    "location": "File > Export > enGrid (.begc)",
    "description": "Export objects as boundaries to enGrid's Blender exchange format (*.begc)",
    "warning": "",
    "wiki_url": "http://engits.eu/wiki",
    "tracker_url": "http://sourceforge.net/apps/mantisbt/engrid",
    "category": "Import-Export"}

'''
Usage Notes:
'''

import bpy
from bpy.props import *
import mathutils, math, struct
from os import remove
import time
from bpy_extras.io_utils import ExportHelper


def do_export(context, props, filepath):
    out = open(filepath, "w")
    N = 0
    for obj in bpy.context.selected_objects:
        if obj.type == 'MESH':
            N = N + 1
        
    out.write('%d\n' % N)
    node_offset = 0    
    
    for obj in bpy.context.selected_objects:
        if obj.type == 'MESH':
            out.write(obj.name)
            out.write('\n')    
    
    for obj in bpy.context.selected_objects:
        if obj.type == 'MESH':
            mesh = obj.data
            mesh.update(calc_tessface=True)
            M = obj.matrix_world
            mesh.transform(M)
            faces = mesh.tessfaces
            nodes = mesh.vertices
            out.write('%d' % len(nodes))
            out.write(' %d\n' % len(faces))
            for n in nodes:
                out.write("%e " % n.co[0])
                out.write("%e " % n.co[1])
                out.write("%e\n" % n.co[2])
            for f in faces:
                out.write("%d" % len(f.vertices))
                for v in f.vertices:
                    out.write(' %d' % (v + node_offset))
                out.write('\n')
            node_offset = node_offset + len(nodes)
            M.invert()
            mesh.transform(M)
    
    out.flush()
    out.close()
    return True


###### EXPORT OPERATOR #######
class Export_engrid(bpy.types.Operator, ExportHelper):
    bl_idname = "export_shape.engrid"
    bl_label = "Export enGrid (.begc)"
    filename_ext = ".begc"
    

    
    @classmethod
    def poll(cls, context):
        return context.active_object.type in {'MESH', 'CURVE', 'SURFACE', 'FONT'}

    def execute(self, context):
        start_time = time.time()
        print('\n_____START_____')
        props = self.properties
        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, self.filename_ext)

        exported = do_export(context, props, filepath)
        
        if exported:
            print('finished export in %s seconds' %((time.time() - start_time)))
            print(filepath)
            
        return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager

        if True:
            # File selector
            wm.fileselect_add(self) # will run self.execute()
            return {'RUNNING_MODAL'}
        elif True:
            # search the enum
            wm.invoke_search_popup(self)
            return {'RUNNING_MODAL'}
        elif False:
            # Redo popup
            return wm.invoke_props_popup(self, event)
        elif False:
            return self.execute(context)


### REGISTER ###

def menu_func(self, context):
    self.layout.operator(Export_engrid.bl_idname, text="enGrid (.begc)")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_export.append(menu_func)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_export.remove(menu_func)
    
if __name__ == "__main__":
    register()
