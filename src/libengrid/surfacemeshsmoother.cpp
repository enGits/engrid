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
#include "surfacemeshsmoother.h"

SurfaceMeshSmoother::SurfaceMeshSmoother()
{
  m_UseEstimatedPlane = false;
  m_Cad = NULL;
  m_Prec = 1e-3;
  m_NumSteps = 5;
  m_UseSimpleCentreScheme = false;
}

void SurfaceMeshSmoother::prepareEstimatedPlane()
{
  m_N0 = vec3_t(0, 0, 0);
  m_X0 = vec3_t(0, 0, 0);
  int count = 0;
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid)) {
      m_N0 += GeometryTools::cellNormal(m_Grid, id_cell);
      m_X0 += cellCentre(m_Grid, id_cell);
      ++count;
    }
  }
  m_N0.normalise();
  m_X0 *= 1.0/count;
  m_UseEstimatedPlane = true;
}

void SurfaceMeshSmoother::prepareCadInterface(CadInterface *cad)
{
  m_Cad = cad;
  m_UseEstimatedPlane = false;
}

vec2_t SurfaceMeshSmoother::transform32(vec3_t x)
{
  x -= m_X0;
  vec3_t x2 = m_M32*x;
  return vec2_t(x2[0], x2[1]);
}

vec3_t SurfaceMeshSmoother::transform23(vec2_t x)
{
  vec3_t x3(x[0], x[1], 0);
  return m_X0 + m_M23*x3;
}

vec3_t SurfaceMeshSmoother::smoothNode(vtkIdType id_node)
{
  // establish 2D coordinate system
  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  if (!m_UseEstimatedPlane) {
    m_X0 = m_Cad->snapNode(id_node, x_old);
    m_N0 = m_Cad->getLastNormal();
  }
  m_N0.normalise();
  m_M32[0] = GeometryTools::orthogonalVector(m_N0);
  m_M32[1] = m_N0.cross(m_M32[0]);
  m_M32[2] = m_N0;
  m_M23 = m_M32.inverse();

  // only smooth simple vertices for now
  //EG_STOPDATE("2015-06-01");
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  if (node_type->GetValue(id_node) != EG_SIMPLE_VERTEX) {
    return x_old;
  }

  // build limiting polygon
  m_Limit.clear();
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node, i);
    EG_GET_CELL(id_cell, m_Grid);
    for (int j = 0; j < num_pts; ++j) {
      vtkIdType id_node1 = pts[num_pts - 1];
      vtkIdType id_node2 = pts[j];
      if (j > 0) {
        id_node1 = pts[j - 1];
      }
      if (id_node1 != id_node && id_node2 != id_node) {
        vec3_t x1, x2;
        m_Grid->GetPoint(id_node1, x1.data());
        m_Grid->GetPoint(id_node2, x2.data());
        m_Limit << QPair<vec2_t,vec2_t>(transform32(x1), transform32(x2));
      }
    }
  }

  vec2_t x(0, 0);
  if (m_UseSimpleCentreScheme) {

    vec2_t v(0, 0);
    for (int i = 0; i < m_Limit.size(); ++i) {
      v += 0.5*(m_Limit[i].first + m_Limit[i].second);
    }
    vec2_t c = v;
    c *= 1.0/m_Limit.size();
    v.normalise();
    if (checkVector(v)) {
      double k1, k2;
      computeLimits(x, v, k1, k2);
      if (k1 > 0) {
        double w = 0.9;
        x = (w*k2 + (1-w)*k1)*v;
      } else {
        x = c;
      }
    }

  } else {

    // compute smallest segment length
    m_MinLength = EG_LARGE_REAL;
    for (int i = 0; i < m_Limit.size(); ++i) {
      vec2_t a = m_Limit[i].first;
      vec2_t b = m_Limit[i].second;
      m_MinLength = min(m_MinLength, (a - b).abs());
    }

    int count = 0;
    vec2_t x1(0,0,0);
    vec2_t x2 = x1;
    do {
      x1 = x2;
      x2 = iteration(x1);
      ++count;
    } while ((x2 - x1).abs() > m_Prec*m_MinLength && count < 10);
    x = x2;

  }

  vec3_t x_smooth = transform23(x);
  if (!m_UseEstimatedPlane) {
    x_smooth = m_Cad->snapNode(id_node, x_smooth);
  }
  return x_smooth;
}

vec2_t SurfaceMeshSmoother::iteration(vec2_t x)
{
  vec2_t v = gradError(x);
  v.normalise();
  if (!checkVector(v)) {
    return x;
  }
  double k1, k2;
  computeLimits(x, v, k1, k2);
  double l1 = k1;
  double l2 = k2;
  while (l2 - l1 > m_Prec) {
    int i = 0;
    double dk = (l2 - l1)/m_NumSteps;
    while (i < m_NumSteps - 2) {
      double j1 = l1 + i*dk;
      double j2 = j1 + dk;
      double err1 = error(x + j1*v);
      double err2 = error(x + j2*v);
      if (err2 > err1) {
        l1 = j1;
        l2 = j2;
        break;
      }
      ++i;
    }
    l1 += i*dk;
    l2  = l1 + dk;
  }
  return x + 0.5*(l1 + l2)*v;
}

double SurfaceMeshSmoother::error(vec2_t x)
{
  if (m_Limit.size() == 0) {
    EG_BUG;
  }

  // min and max area
  double A_min =  EG_LARGE_REAL;
  double A_max = -EG_LARGE_REAL;
  for (int i = 0; i < m_Limit.size(); ++i) {
    vec2_t a = m_Limit[i].first - x;
    vec2_t b = m_Limit[i].second - x;
    double A = a[0]*b[1] - a[1]*b[0];
    A_min = min(A_min, A);
    A_max = max(A_max, A);
  }
  double err = A_max - A_min;
  return err;
  /*
  //A0 /= m_Limit.count();

  // compute error
  double err = 0;
  for (int i = 0; i < m_Limit.size(); ++i) {
    vec2_t a = m_Limit[i].first - x;
    vec2_t b = m_Limit[i].second - x;
    double A = a[0]*b[1] - a[1]*b[0];
    err += sqr(A - A0);
  }

  return err;
  */
}

vec2_t SurfaceMeshSmoother::gradError(vec2_t x)
{
  double d = 1e-3*m_MinLength;
  vec2_t grad;
  grad[0] = (error(x + vec2_t(0.5*d, 0)) - error(x - vec2_t(0.5*d, 0)))/d;
  grad[1] = (error(x + vec2_t(0, 0.5*d)) - error(x - vec2_t(0, 0.5*d)))/d;
  return grad;
}

void SurfaceMeshSmoother::computeLimits(vec2_t x, vec2_t v, double &k1, double &k2)
{
  int count = 0;
  k1 = EG_LARGE_REAL;
  k2 = -k1;
  for (int i = 0; i < m_Limit.size(); ++i) {
    vec2_t a = m_Limit[i].first;
    vec2_t b = m_Limit[i].second;
    double k, l;
    GeometryTools::intersection(k, l, x, v, a, (b - a));
    if (l >= 0 && l <= 1) {
      k1 = min(k1, k);
      k2 = max(k2, k);
      ++count;
    }
  }
  if (count < 2) {
    //EG_BUG;
  }
}

void SurfaceMeshSmoother::operate()
{
  EG_BUG;
}
