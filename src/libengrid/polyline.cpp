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
#include "polyline.h"

void PolyLine::getStencil(double l, int &i1, int &i2, double &w1, double &w2)
{
  l *= m_Length;
  double L = 0;
  w2 = 1;
  i1 = m_Points.size() - 2;
  i2 = m_Points.size() - 1;
  for (int i = 1; i < m_Points.size(); ++i) {
    double dL = (m_Points[i] - m_Points[i-1]).abs();
    if (l >= L && l <= L + dL) {
      i1 = i - 1;
      i2 = i;
      w2 = (l - L)/dL;
      break;
    }
    L += dL;
  }
  w1 = 1 - w2;
}

vec3_t PolyLine::position(double l)
{
  int i1, i2;
  double w1, w2;
  getStencil(l, i1, i2, w1, w2);
  return w1*m_Points[i1] + w2*m_Points[i2];
}

vec3_t PolyLine::normal(double l)
{
  int i1, i2;
  double w1, w2;
  getStencil(l, i1, i2, w1, w2);
  vec3_t n = m_Points[i2] - m_Points[i1];
  n.normalise();
  return n;
}

vec3_t PolyLine::intersection(vec3_t x, vec3_t n)
{
  double min_err = EG_LARGE_REAL;
  vec3_t x_best = x;
  for (int i = 1; i < m_Points.size(); ++i) {
    double k = GeometryTools::intersection(m_Points[i-1], m_Points[i] - m_Points[i-1], x, n);
    double err = fabs(k - 0.5);
    if (err < min_err) {
      min_err = err;
      x_best = (1-k)*m_Points[i-1] + k*m_Points[i];
    }
  }
  return x_best;
}
