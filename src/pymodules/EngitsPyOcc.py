#!/usr/bin/python
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

from OCC.Utils import Common, Construct
from OCC.StlAPI import *
from OCC.DataExchange.STEP import *
from OCC.DataExchange.IGES import *

import sys
import math

def angle3D(ux, uy, uz, vx, vy, vz):
  ptot2 = (ux**2 + uy**2 + uz**2)*(vx**2 + vy**2 + vz**2)
  if ptot2 > 0:
    arg = (ux*vx + uy*vy + uz*vz)/math.sqrt(ptot2)
    if arg > 1.0:
      arg =  1.0
    if arg < -1.0: 
      arg = -1.0
    return math.degrees(math.acos(arg))
  return 0.0


def partial_cylinder(x1, y1, z1, x2, y2, z2, radius, angle):
  L = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
  direction = Common.gp_Dir((x2 - x1)/L, (y2 - y1)/L, (z2 - z1)/L)
  axis = Common.gp_Ax2(Common.gp_Pnt(x1, y1, z1), direction)
  if angle < 360:
    angle = angle*math.pi/180
    shape = BRepPrimAPI.BRepPrimAPI_MakeCylinder(axis, radius, L, angle).Shape()
  else:
    shape = BRepPrimAPI.BRepPrimAPI_MakeCylinder(axis, radius, L).Shape()
  return shape


def cylinder(x1, y1, z1, x2, y2, z2, radius):
  return partial_cylinder(x1, y1, z1, x2, y2, z2, radius, 360)


def box(x1, y1, z1, x2, y2, z2):
  return BRepPrimAPI.BRepPrimAPI_MakeBox(Common.gp_Pnt(x1, y1, z1), Common.gp_Pnt(x2, y2, z2)).Shape()


def rotate(shape, x1, y1, z1, x2, y2, z2, angle):
  L = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
  direction = Common.gp_Dir((x2 - x1)/L, (y2 - y1)/L, (z2 - z1)/L)
  axis = Common.gp_Ax1(Common.gp_Pnt(x1, y1, z1), direction)
  return Construct.rotate(shape, axis, angle)


def rotate_x(shape, x, y, z, angle):
  return rotate(shape, x, y, z, x + 1, y, z, angle)


def rotate_x0(shape, angle):
  return rotate(shape, 0, 0, 0, 1, 0, 0, angle)


def rotate_y(shape, x, y, z, angle):
  return rotate(shape, x, y, z, x, y + 1, z, angle)


def rotate_y0(shape, angle):
  return rotate(shape, 0, 0, 0, 0, 1, 0, angle)


def rotate_z(shape, x, y, z, angle):
  return rotate(shape, x, y, z, x, y, z + 1, angle)


def rotate_z0(shape, angle):
  return rotate(shape, 0, 0, 0, 0, 0, 1, angle)


def cut(shape1, shape2):
  return Construct.boolean_cut(shape1, shape2)


def join(shape1, shape2):
  return Construct.boolean_fuse(shape1, shape2)


def intersect(shape1, shape2):
  cut_shape = cut(shape1, shape2)
  return cut(shape1, cut_shape)


def move(shape, dx, dy, dz):
  return Construct.translate_topods_from_vector(shape, Common.gp_Vec(dx, dy, dz))


def move_x(shape, dx):
  return Construct.translate_topods_from_vector(shape, Common.gp_Vec(dx, 0, 0))


def move_y(shape, dy):
  return Construct.translate_topods_from_vector(shape, Common.gp_Vec(0, dy, 0))


def move_z(shape, dz):
  return Construct.translate_topods_from_vector(shape, Common.gp_Vec(0, 0, dz))


def cone(x1, y1, z1, x2, y2, z2, R1, R2):
  L = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
  angle = angle3D(0, 0, 1, x2-x1, y2-y1, z2-z1)
  shape = BRepPrimAPI.BRepPrimAPI_MakeCone(R1, R2, L).Shape()
  l = math.sqrt((x2-x1)**2 + (y2-y1)**2)
  if l/L > 1e-4:
    shape = rotate(shape, 0, 0, 0, y1-y2, x2-x1, 0, angle)
  return move(shape, x1, y1, z1)


def sphere(x, y, z, R):
  shape = BRepPrimAPI.BRepPrimAPI_MakeSphere(R).Shape()
  return move(shape, x, y, z)


def nose(x1, y1, z1, x2, y2, z2, R, r):
  nx = x2 - x1
  ny = y2 - y1
  nz = z2 - z1
  H = math.sqrt(nx*nx + ny*ny + nz*nz)
  nx /= H
  ny /= H
  nz /= H
  xe = R;
  ye = r - H
  scal = 1
  a1 = 0
  a2 = 0.5*math.pi
  xt = 0
  yt = 0
  count = 0
  while math.fabs(scal) > 1e-6 and count < 1000:
    a = 0.5*(a1 + a2)
    xt = r*math.cos(a)
    yt = r*math.sin(a)
    vx = xt - xe
    vy = yt - ye
    L = math.sqrt(vx*vx + vy*vy)
    scal = (vx*xt + vy*yt)/(L*r)
    if scal > 0:
      a2 = a
    else:
      a1 = a
    count += 1
  xs = x2 - r*nx
  ys = y2 - r*ny
  zs = z2 - r*nz
  xc = xs + yt*nx
  yc = ys + yt*ny
  zc = zs + yt*nz
    
  cone_shape = cone(x1, y1, z1, xc, yc, zc, R, xt)
  sphere_shape = sphere(xs, ys, zs, r)
   
  return join(cone_shape, sphere_shape)


def scale(shape, x, y, z, factor):
  return Construct.scale(shape, Common.gp_Pnt(x,y,z), factor)


def write_stl(shape, file_name, precision):
  writer = StlAPI_Writer()
  writer.SetRelativeMode(True)
  writer.SetCoefficient(precision)
  writer.Write(shape, file_name)



