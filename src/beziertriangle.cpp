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
#include "beziertriangle.h"
#include "engrid.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include "vtkUnstructuredGridWriter.h"

#include <vtkCellLocator.h>

BezierTriangle::BezierTriangle()
{
}

BezierTriangle::BezierTriangle(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110)
{
  setControlPoints( X_200,  X_020,  X_002,  X_011,  X_101,  X_110);
}

void BezierTriangle::setControlPoints(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110)
{
  m_X_200 = X_200;
  m_X_020 = X_020;
  m_X_002 = X_002;
  m_X_011 = X_011;
  m_X_101 = X_101;
  m_X_110 = X_110;
}

void BezierTriangle::getControlPoints(vec3_t& X_200, vec3_t& X_020, vec3_t& X_002, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110)
{
  X_200 = m_X_200;
  X_020 = m_X_020;
  X_002 = m_X_002;
  X_011 = m_X_011;
  X_101 = m_X_101;
  X_110 = m_X_110;
}

void BezierTriangle::writeBezierSurface()
{
  //qDebug()<<"writeBezierSurface called";
  int N=10;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  
  //qDebug()<<"N_cells="<<N_cells;
  //qDebug()<<"N_points="<<N_points;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, 2*N_cells, 2*N_points);
  
  vtkIdType offset = 0;
  offset += addBezierSurface(this, bezier, offset, N);
  
//   BezierTriangle B(m_X_200, m_X_020, m_X_002, m_X_011-vec3_t(0,0,1), m_X_101-vec3_t(0,0,1), m_X_110-vec3_t(0,0,1));
//   offset += B.addBezierSurface(bezier, offset, N);
  
  //qDebug()<<"offset="<<offset;
  
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu2);
  vtu2->SetFileName("bezier.vtu");
  vtu2->SetDataModeToBinary();
//   vtu2->SetDataModeToAscii();
  vtu2->SetInput(bezier);
  vtu2->Write();
  
}

vec3_t BezierTriangle::QuadraticBezierTriangle(double u, double v, double w)
{
  double total = u + v + w;
  u=u/total;
  v=v/total;
  w=w/total;
  return pow(u,2)*m_X_200 + pow(v,2)*m_X_020 + pow(w,2)*m_X_002 + 2*u*v*m_X_110 + 2*v*w*m_X_011 + 2*w*u*m_X_101;
}

vec3_t BezierTriangle::QuadraticBezierTriangle(vec2_t M)
{
  vec3_t bary_coords = getBarycentricCoordinates(M[0],M[1]);
  double u,v,w;
  u=bary_coords[0];
  v=bary_coords[1];
  w=bary_coords[2];
  return QuadraticBezierTriangle(u, v, w);
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle(vec3_t g_M)
{
  int N=10;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, 2*N_cells, 2*N_points);
  
  vtkIdType offset = 0;
  offset += addBezierSurface(this, bezier, offset, N);
  
  vtkIdType cellId;
  int subId;
  double dist2;
  vtkCellLocator* locator=vtkCellLocator::New();
  locator->SetDataSet(bezier);
  locator->BuildLocator();
  vec3_t g_P;
  locator->FindClosestPoint(g_M.data(),g_P.data(),cellId,subId,dist2);
  locator->Delete();
  return g_P;
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle2(vec3_t g_M)
{
  double L = (m_X_200-m_X_020).abs();
  L = min(L,(m_X_020-m_X_002).abs());
  L = min(L,(m_X_002-m_X_200).abs());
  
  double maxerr = L/100.;
  vec3_t xi;
  vec3_t ri;
  double d;
  projectOnTriangle(g_M, xi, ri, d);
  vec3_t g_A = xi;
  vec2_t t_A = vec2_t(ri[0],ri[1]);
  vec3_t g_B = QuadraticBezierTriangle(t_A);
  projectOnTriangle(g_B, xi, ri, d);
  vec3_t g_C = xi;
  vec3_t err_vector = g_A - g_C;
  while(err_vector.abs()>maxerr) {
    g_A = g_A + err_vector;
    projectOnTriangle(g_A, xi, ri, d);
    t_A = vec2_t(ri[0],ri[1]);
    g_B = QuadraticBezierTriangle(t_A);
    projectOnTriangle(g_B, xi, ri, d);
    vec3_t g_C = xi;
    err_vector = g_A - g_C;
  }
  return g_B;
}

bool BezierTriangle::projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d)
{
  vec3_t T_a = m_X_200;
  vec3_t T_b = m_X_020;
  vec3_t T_c = m_X_002;
  vec3_t T_g1 = T_b-T_a;
  vec3_t T_g2 = T_c-T_a;
  vec3_t T_g3 = T_g1.cross(T_g2);
  
  xi = vec3_t(1e99,1e99,1e99);
  double scal = (xp - T_a)*T_g3;
  vec3_t x1, x2;
  if (scal > 0) {
    x1 = xp + T_g3;
    x2 = xp - scal*T_g3 - T_g3;
  } else {
    x1 = xp - T_g3;
    x2 = xp - scal*T_g3 + T_g3;
  }
  d = 1e99;
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T_a, T_b, T_c, x1, x2, xi, ri);
  if (intersects_face) {
    vec3_t dx = xp - T_a;
    d = fabs(dx*T_g3);
  } else {
    double kab = GeometryTools::intersection(T_a, T_b - T_a, xp, T_b - T_a);
    double kac = GeometryTools::intersection(T_a, T_c - T_a, xp, T_c - T_a);
    double kbc = GeometryTools::intersection(T_b, T_c - T_b, xp, T_c - T_b);
    double dab = (T_a + kab*(T_b-T_a) - xp).abs();
    double dac = (T_a + kac*(T_c-T_a) - xp).abs();
    double dbc = (T_b + kbc*(T_c-T_b) - xp).abs();
    bool set = false;
    if ((kab >= 0) && (kab <= 1)) {
      if (dab < d) {
        xi = T_a + kab*(T_b-T_a);
        d = dab;
        set = true;
      }
    }
    if ((kac >= 0) && (kac <= 1)) {
      if (dac < d) {
        xi = T_a + kac*(T_c-T_a);
        d = dac;
        set = true;
      }
    }
    if ((kbc >= 0) && (kbc <= 1)) {
      if (dbc < d) {
        xi = T_b + kbc*(T_c-T_b);
        d = dbc;
        set = true;
      }
    }
    double da = (T_a - xp).abs();
    double db = (T_b - xp).abs();
    double dc = (T_c - xp).abs();
    if (da < d) {
      xi = T_a;
      d = da;
      set = true;
    }
    if (db < d) {
      xi = T_b;
      d = db;
    }
    if (dc < d) {
      xi = T_c;
      d = dc;
      set = true;
    }
    if (!set) {
      EG_BUG;
    }
  }
  if (xi[0] > 1e98) {
    EG_BUG;
  }
  return intersects_face;
}
