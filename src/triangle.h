//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vtkIdList.h"
#include "vtkUnstructuredGrid.h"
#include "math/mathvector.h"
#include "math/smallsquarematrix.h"

class Triangle{
public:
    vtkIdType id_a, id_b, id_c;
    vec3_t a, b, c;
    vec3_t g1, g2, g3;
    mat3_t G, GI;
    double A;
    double smallest_length;
  
public:
  Triangle();
  Triangle(vec3_t a_a, vec3_t a_b, vec3_t a_c);
  Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_a, vtkIdType a_id_b, vtkIdType a_id_c);
  Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_cell);
  void setupTriangle();
  
public:
  bool projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d);
  vec3_t local3DToGlobal3D(vec3_t l_M);
  vec3_t global3DToLocal3D(vec3_t g_M);
  vec3_t local2DToGlobal3D(vec2_t l_M);
  vec2_t global3DToLocal2D(vec3_t g_M);
};

#endif
