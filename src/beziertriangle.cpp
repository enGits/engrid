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

#include "vtkUnstructuredGridWriter.h"

int idx_func(int N, int i, int j)
{
  int offset = -i*(i-2*N-1)/2;
  return offset+j;
}

vec3_t getBarycentricCoordinates(double x, double y)
{
  // initialize
  double t1=0;
  double t2=0;
  double t3=0;
  
/*  if(x==0) {
    t3=y;
    t1=1-y;
    t2=0;
  }
  else if(y==0) {
    t2=x;
    t1=1-x;
    t3=0;
  }
  else if((x+y)==1) {
  
  }
  else {
  }
  
  double k1,k2;
  if(!intersection (k1, k2, t_A, t_M-t_A, t_B, t_C-t_B)) EG_BUG;
  vec2_t t_I1 = t_A+k1*(t_M-t_A);
  vec3_t g_nI1 = (1-k2)*g_nB + k2*g_nC;
  vec2_t pm1_M(1.0/k1,0);
  
  // normalize
  double total = t1+t2+t3;
  t1=t1/total;
  t2=t2/total;
  t3=t3/total;*/
  
  t2 = x;
  t3 = y;
  t1 = 1-t2-t3;
  
  // return value
  vec3_t bary_coords(t1,t2,t3);
  return bary_coords;
}

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
  offset += addBezierSurface(bezier, offset, N);
  
  BezierTriangle B(m_X_200, m_X_020, m_X_002, m_X_011-vec3_t(0,0,1), m_X_101-vec3_t(0,0,1), m_X_110-vec3_t(0,0,1));
  offset += B.addBezierSurface(bezier, offset, N);
  
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

vtkIdType BezierTriangle::addBezierSurface(vtkUnstructuredGrid* bezier, int offset, int N)
{
  vtkIdType node_count = 0;
  for(int i=0;i<N;i++) {
    for(int j=0;j<N-i;j++) {
      double x = i/(double)(N-1);
      double y = j/(double)(N-1);
      vec3_t bary_coords = getBarycentricCoordinates(x,y);
      double u,v,w;
      u=bary_coords[0];
      v=bary_coords[1];
      w=bary_coords[2];
      vec3_t M = QuadraticBezierTriangle(u, v, w);
      bezier->GetPoints()->SetPoint(offset + node_count, M.data());node_count++;
    }
  }
  
//   qDebug()<<"node_count="<<node_count;
  
  int cell_count = 0;
  for(int i=0;i<N-1;i++) {
    for(int j=0;j<N-1-i;j++) {
      
      //qDebug()<<"(i,j)="<<i<<j;
      
      vtkIdType pts_triangle1[3];
      pts_triangle1[0]=offset + idx_func(N, i  ,j  );
      pts_triangle1[1]=offset + idx_func(N, i+1,j  );
      pts_triangle1[2]=offset + idx_func(N, i  ,j+1);
      bezier->InsertNextCell(VTK_TRIANGLE,3,pts_triangle1);cell_count++;
      
      if(i+j<N-2) {
        //qDebug()<<"BEEP";
        vtkIdType pts_triangle2[3];
        pts_triangle2[0]=offset + idx_func(N, i+1,j  );
        pts_triangle2[1]=offset + idx_func(N, i+1,j+1);
        pts_triangle2[2]=offset + idx_func(N, i  ,j+1);
        bezier->InsertNextCell(VTK_TRIANGLE,3,pts_triangle2);cell_count++;
      }
      
    }
  }
  
  //qDebug()<<"cell_count="<<cell_count;
  return node_count;
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

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle(double u, double v, double w)
{
/*  vec3_t B = QuadraticBezierTriangle();
  B*/
}
