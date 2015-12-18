// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "curve.h"

mat3_t Curve::computeBase(double l, Curve* curve)
{
  vec3_t g3 = normal(l);
  g3.normalise();
  vec3_t g1 = curve->position(l) - position(l);
  g1.normalise();
  vec3_t g2 = g3.cross(g1);
  g2.normalise();
  mat3_t A;
  clinit(A) = g1, g2, g3;
  return A.transp();
}

mat3_t Curve::computeOrthoBase(double l, Curve *curve)
{
  vec3_t g3 = normal(l);
  g3.normalise();
  vec3_t x = position(l);
  vec3_t x_curve = curve->intersection(x, g3);
  vec3_t g1 = x_curve - x;
  g1.normalise();
  vec3_t g2 = g3.cross(g1);
  g2.normalise();
  mat3_t A;
  clinit(A) = g1, g2, g3;
  return A.transp();
}
