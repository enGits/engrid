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
    "name": "Import from enGrid (.begc)",
    "author": "Oliver Gloth",
    "version": (0, 1),
    "blender": (2, 5, 9),
    "api": 36079,
    "location": "File > Import > enGrid (.begc)",
    "description": "Import objects from enGrid's Blender exchange format (*.begc)",
    "warning": "",
    "wiki_url": "http://engits.eu/wiki",
    "tracker_url": "http://sourceforge.net/apps/mantisbt/engrid",
    "category": "Import-Export"}

import bpy
import time
from bpy.props import *
from mathutils import *
from os import remove
from bpy_extras.io_utils import *



def do_import(context, props, filepath):
    in_file = open(filepath, "r")
    line = in_file.readline()
    Nobjects = int(line)
    object_names = []
    for i in range(0, Nobjects):
        line = in_file.readline()
        object_names.append(line.strip())
    
    global_verts = []
    offset = 0
    for i_object in range(0, Nobjects):
        line = in_file.readline()
        words = line.split()
        Nverts = int(words[0])
        Nfaces = int(words[1])
        
        local_verts = []
        for i_vert in range(0, Nverts):
            line = in_file.readline()
            words = line.split()
            x = float(words[0])
            y = float(words[1])
            z = float(words[2])
            local_verts.append( Vector((x,y,z)) )
            global_verts.append( Vector((x,y,z)) )
        
        faces = []
        print ("Nfaces=", Nfaces)
        for i_face in range(0, Nfaces):
            line = in_file.readline()
            words = line.split()
            if len(words) < 3:
                return
            Nverts_in_face = int(words[0])
            if len(words) != 1 + Nverts_in_face:
                return
            face_verts = []
            for i_face_vert in range(0, Nverts_in_face):
                idx = int(words[i_face_vert + 1]) - offset
                face_verts.append(idx)
            faces.append(face_verts)
        
        mesh = bpy.data.meshes.new(object_names[i_object])
        mesh.from_pydata(local_verts, [], faces)
        mesh.update()
        from bpy_extras import object_utils
        object_utils.object_data_add(context, mesh, operator=None)
        #BPyAddMesh.add_mesh_simple(object_names[i_object], local_verts, [], faces)

        offset += Nverts

    in_file.close()
    return True

    
    
###### IMPORT OPERATOR #######
class Import_engrid(bpy.types.Operator, ImportHelper):
    bl_idname = "import_shape.engrid"
    bl_label = "Import enGrid (.begc)"
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
        do_import(context, props, filepath)
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
    self.layout.operator(Import_engrid.bl_idname, text="enGrid (.begc)")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_func)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_func)
    
if __name__ == "__main__":
    register()

