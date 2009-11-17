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
#include "triangle.h"

#include <vtkUnstructuredGrid.h>

class BezierTriangle : public Triangle, public EgVtkObject
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
  void writeBezierSurface(QString filename, int N);
  
  vec3_t QuadraticBezierTriangle(double u, double v, double w);
  vec3_t QuadraticBezierTriangle(vec2_t M);
  vec3_t QuadraticBezierTriangle_g(vec3_t g_M);
  
  vec3_t projectOnQuadraticBezierTriangle(vec3_t g_M);
  vec3_t projectOnQuadraticBezierTriangle2(vec3_t g_M);
  vec3_t projectOnQuadraticBezierTriangle3(vec3_t g_M, int output=0);
  vec3_t projectOnQuadraticBezierTriangle4(vec3_t g_M);
  vec3_t projectOnQuadraticBezierTriangle5(vec3_t g_M);
  
// stuff used for projections on the Bezier surface
private:
  void setupFunctionVariables();
  
  vec3_t m_l_X_200;
  vec3_t m_l_X_020;
  vec3_t m_l_X_002;
  vec3_t m_l_X_011;
  vec3_t m_l_X_101;
  vec3_t m_l_X_110;
  
  vec3_t m_coeff_x2;
  vec3_t m_coeff_y2;
  vec3_t m_coeff_xy;
  vec3_t m_coeff_x;
  vec3_t m_coeff_y;
public:
  vec2_t fixedPointFunction(vec2_t t_inputPoint, double x, double y);
  vec2_t fixedPointFunction(vec2_t t_inputPoint, vec2_t A);
  mat2_t jacobiMatrix(double x, double y);
  mat2_t jacobiMatrix_numeric(vec2_t t_inputPoint, double x, double y, double dx, double dy);
  
//   mat3_t jacobiMatrix_no_projection(double x, double y);
  vec3_t surfaceNormal(vec2_t t_M, int output);
  double z_func(vec2_t t_M);
  double z_func(double x, double y);
};

#endif
