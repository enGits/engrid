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

BezierTriangle::BezierTriangle() : Triangle(), EgVtkObject()
{
}

BezierTriangle::BezierTriangle(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110)  : Triangle(X_200,X_020,X_002), EgVtkObject()
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
  setupFunctionVariables();
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
  int N=100;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  
  //qDebug()<<"N_cells="<<N_cells;
  //qDebug()<<"N_points="<<N_points;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, N_cells, N_points);
  
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
  int N=100;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, N_cells, N_points);
  
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
  vec3_t err_vector = g_M - g_C;
  int Nloop=0;
  while(err_vector.abs()>maxerr && Nloop<200) {
//   for(int i=0;i<10;i++) {
    g_A = g_A + err_vector;
    projectOnTriangle(g_A, xi, ri, d);
    t_A = vec2_t(ri[0],ri[1]);
    g_B = QuadraticBezierTriangle(t_A);
    projectOnTriangle(g_B, xi, ri, d);
    vec3_t g_C = xi;
    err_vector = g_M - g_C;
    Nloop++;
  }
  return g_B;
}

void BezierTriangle::setupFunctionVariables() {
  m_t_X_200 = globalToLocal(m_X_200);
  m_t_X_020 = globalToLocal(m_X_020);
  m_t_X_002 = globalToLocal(m_X_002);
  m_t_X_011 = globalToLocal(m_X_011);
  m_t_X_101 = globalToLocal(m_X_101);
  m_t_X_110 = globalToLocal(m_X_110);
  
  m_coeff_x2 = m_t_X_020 - 2*m_t_X_110;
  m_coeff_y2 = m_t_X_002 - 2*m_t_X_101;
  m_coeff_xy = -2*m_t_X_110 + 2*m_t_X_011 - 2*m_t_X_101;
  m_coeff_x = 2*m_t_X_110;
  m_coeff_y = 2*m_t_X_101;
}

vec2_t BezierTriangle::fixedPointFunction(vec2_t t_inputPoint, double x, double y)
{
  vec2_t F;
  F[0] = pow(x,2)*m_coeff_x2[0] + pow(y,2)*m_coeff_y2[0] + x*y*m_coeff_xy[0] + x*m_coeff_x[0] + y*m_coeff_y[0] - t_inputPoint[0];
  F[1] = pow(x,2)*m_coeff_x2[1] + pow(y,2)*m_coeff_y2[1] + x*y*m_coeff_xy[1] + x*m_coeff_x[1] + y*m_coeff_y[1] - t_inputPoint[1];
  return F;
}

mat2_t BezierTriangle::jacobiMatrix(double x, double y)
{
  mat2_t J;
  J[0][0] = 2*x*m_coeff_x2[0] + y*m_coeff_xy[0] + m_coeff_x[0];
  J[1][0] = 2*x*m_coeff_x2[1] + y*m_coeff_xy[1] + m_coeff_x[1];
  J[0][1] = 2*y*m_coeff_y2[0] + x*m_coeff_xy[0] + m_coeff_y[0];
  J[1][1] = 2*y*m_coeff_y2[1] + x*m_coeff_xy[1] + m_coeff_y[1];
  return J;
}
