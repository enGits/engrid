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
  int N=2;
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
  
/*//   EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  EG_VTKSP(vtkUnstructuredGridWriter,vtu);
  vtu->SetFileName("bezier.vtk");
//   vtu->SetDataModeToBinary();
//   vtu->SetDataModeToAscii();
  vtu->SetInput(bezier);
  vtu->Write();*/
  
  EG_VTKSP(vtkUnstructuredGridWriter,vtu1);
  vtu1->SetFileName("bezier.vtk");
  vtu1->SetInput(bezier);
  vtu1->Write();
  
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
