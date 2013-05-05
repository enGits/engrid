# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Import Selig airfoil (.dat)",
    "author": "Oliver Gloth (enGits GmbH)",
    "version": (1, 0),
    "blender": (2, 6, 3),
    "location": "File > Import > Selig airfoil (.dat)",
    "description": "Import edges from a Selig airfoil file (*.dat)",
    "warning": "",
    "wiki_url": "https://github.com/enGits/engrid/wiki",
    "tracker_url": "https://github.com/enGits/engrid/issues",
    "category": "Import-Export"}

import bpy
import time
from bpy.props import *
from mathutils import *
from os import remove
from bpy_extras.io_utils import *



def do_import1(context, props, filepath):
    in_file = open(filepath, "r")
    object_name = in_file.readline()

    line = in_file.readline()
    x = []
    y = []
    z = []
    while line:
        words = line.split()
        x.append(float(words[0]))
        y.append(float(words[1]))
        z.append(float(0))
        line = in_file.readline()        

    verts = []
    for i in range(0, len(x)):
        verts.append(Vector((x[i],y[i],z[i])))
    verts.append(Vector((x[0],y[0],z[0])))
    
    edges = []

    for i in range(0, len(verts) - 1):
        edge_verts = []
        edge_verts.append(i)
        edge_verts.append(i + 1)
        edges.append(edge_verts)
    
    mesh = bpy.data.meshes.new(object_name)
    mesh.from_pydata(verts, edges, [])
    mesh.update()
    from bpy_extras import object_utils
    object_utils.object_data_add(context, mesh, operator=None)

    in_file.close()
    return True
    
def save_read(in_file):
    line = in_file.readline()
    while len(line) < 2:
        line = in_file.readline()
    return line

def do_import2(context, props, filepath):
    in_file = open(filepath, "r")
    object_name = save_read(in_file)

    line = save_read(in_file)
    words = line.split()
    N = []
    N.append(int(float(words[0])))
    N.append(int(float(words[1])))
        
    x = []
    y = []
    z = []
    
    verts = []
    for i in range(2):
        for j in range(N[i]):
            line = save_read(in_file)
            words = line.split()
            verts.append(Vector((float(words[0]), float(words[1]), float(0))))
            
    edges = []
    for i in range(2):
        for j in range(N[i] - 1):
            edge_verts = []
            edge_verts.append(i*N[0] + j)
            edge_verts.append(i*N[0] + j + 1)
            edges.append(edge_verts)
        
    mesh = bpy.data.meshes.new(object_name)
    mesh.from_pydata(verts, edges, [])
    mesh.update()
    from bpy_extras import object_utils
    object_utils.object_data_add(context, mesh, operator=None)

    in_file.close()
    return True

   
    
###### IMPORT OPERATOR #######
class Import_selig(bpy.types.Operator, ImportHelper):
    bl_idname = "import_shape.selig"
    bl_label = "Import Selig airfoil (.dat)"
    filename_ext = ".dat"

    loop = BoolProperty(name="node-loop", description="Are all nodes in a loop?", default=True)

    @classmethod
    def poll(cls, context):
        return context.active_object.type in {'MESH', 'CURVE', 'SURFACE', 'FONT'}

    def execute(self, context):
        start_time = time.time()
        print('\n_____START_____')
        props = self.properties
        filepath = self.filepath
        filepath = bpy.path.ensure_ext(filepath, self.filename_ext)
        if self.loop:
            do_import1(context, props, filepath)
        else:
            do_import2(context, props, filepath)
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
    self.layout.operator(Import_selig.bl_idname, text="Selig airfoil (.dat)")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_func)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_func)
    
if __name__ == "__main__":
    register()

