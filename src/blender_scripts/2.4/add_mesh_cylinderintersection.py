#!BPY
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
"""
Name: 'Cylinder Intersection'
Blender: 249
Group: 'AddMesh'
"""
import BPyAddMesh
import Blender
try: from math import cos, sin, pi, sqrt
except: math = None

def add_cylinters(PREF_MAJOR_RAD, PREF_MINOR_RAD, PREF_OFFSET, PREF_SEG):
	Vector = Blender.Mathutils.Vector
	verts = []
	faces = []
	i1 = 0
	tot_verts = PREF_SEG
	verts_tmp = []
	
	for i in xrange(PREF_SEG):
		angle = 2*pi*i/PREF_SEG
		
		y = (cos(angle)*PREF_MINOR_RAD) + PREF_OFFSET
		z = (sin(angle)*PREF_MINOR_RAD)
		x = -sqrt(PREF_MAJOR_RAD**2 - y**2)
		verts.append(Vector(x, y, z))
		if i+1 == PREF_SEG:
			i2 = 0
		else:
			i2 = i1 + 1
		
		faces.append((i1,i2))
		i1 += 1

	return verts, faces

def main():
	Draw = Blender.Draw
	PREF_MAJOR_RAD = Draw.Create(1.0)
	PREF_MINOR_RAD = Draw.Create(0.1)
	PREF_OFFSET = Draw.Create(0.0)
	PREF_SEG = Draw.Create(36)

	if not Draw.PupBlock('Add Cylinder Intersection', [\
	('Major Radius:', PREF_MAJOR_RAD,  0.001, 10000, 'Radius for the bigger cylinder'),\
	('Minor Radius:', PREF_MINOR_RAD,  0.001, 10000, 'Radius for the smaller cylinder'),\
	('Axis Offset:', PREF_OFFSET,  0.000, 10000, 'Axix offset between smaller and bigger cylinder'),\
	('Segments:', PREF_SEG,  3, 256, 'Number of segments for the intersection (circle)'),\
	]):
		return
	
	verts, faces = add_cylinters(PREF_MAJOR_RAD.val, PREF_MINOR_RAD.val, PREF_OFFSET.val, PREF_SEG.val)
	
	BPyAddMesh.add_mesh_simple('CylInters', verts, [], faces)

if cos and sin and pi:
    main()
else:
    Blender.Draw.PupMenu("Error%t|This script requires a full python installation")

