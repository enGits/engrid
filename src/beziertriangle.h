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
#ifndef BEZIERTRIANGLE_H
#define BEZIERTRIANGLE_H

#include "math/mathvector.h"
#include "math/linsolve.h"
#include "math/smallsquarematrix.h"

#include "egvtkobject.h"

#include <vtkUnstructuredGrid.h>

class BezierTriangle : public EgVtkObject
{
public:
  vec3_t m_X_200;
  vec3_t m_X_020;
  vec3_t m_X_002;
  vec3_t m_X_011;
  vec3_t m_X_101;
  vec3_t m_X_110;
  
public:
  BezierTriangle();
  BezierTriangle(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  
  void setControlPoints(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  void getControlPoints(vec3_t& X_200, vec3_t& X_020, vec3_t& X_002, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
  vec2_t BaryCentric2xy(vec3_t barycoords);
  vec3_t xy2BaryCentric(vec2_t xycoords);
  vec3_t Bezier(vec2_t xycoords);
  vec3_t Bezier(vec3_t barycoords);
  vec3_t Projection(vec2_t xycoords);
  vec3_t Projection(vec3_t barycoords);
  mat2_t JacobiMatrix();
  vec2_t FixedPointFunction(vec2_t xycoords);
  void writeBezierSurface();
//   vtkIdType addBezierSurface(vtkUnstructuredGrid* bezier, int offset, int N);
  vec3_t QuadraticBezierTriangle(double u, double v, double w);
  vec3_t QuadraticBezierTriangle(vec2_t M);
  vec3_t projectOnQuadraticBezierTriangle(vec3_t M);
  vec3_t projectOnQuadraticBezierTriangle2(vec3_t M);
  
  bool projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d);
  
};

#endif
