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

void BezierTriangle::writeBezierSurface(QString filename, int N)
{
  //qDebug()<<"writeBezierSurface called";
//   int N=10;
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
  vtu2->SetFileName(qPrintable(filename));
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

vec3_t BezierTriangle::QuadraticBezierTriangle_g(vec3_t g_M)
{
  vec2_t t_M = global3DToLocal2D(g_M);
  return QuadraticBezierTriangle(t_M);
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
  projectOnTriangle(g_M, xi, ri, d, false);
  vec3_t g_A = xi;
  vec2_t t_A = vec2_t(ri[0],ri[1]);
  vec3_t g_B = QuadraticBezierTriangle(t_A);
  projectOnTriangle(g_B, xi, ri, d, false);
  vec3_t g_C = xi;
  vec3_t err_vector = g_M - g_C;
  int Nloop=0;
/*  while(err_vector.abs()>maxerr && Nloop<200) {
//   for(int i=0;i<10;i++) {
    g_A = g_A + err_vector;
    projectOnTriangle(g_A, xi, ri, d, false);
    t_A = vec2_t(ri[0],ri[1]);
    g_B = QuadraticBezierTriangle(t_A);
    projectOnTriangle(g_B, xi, ri, d, false);
    vec3_t g_C = xi;
    err_vector = g_M - g_C;
    Nloop++;
  }*/
  return g_B;
}

void BezierTriangle::setupFunctionVariables() {
  m_l_X_200 = global3DToLocal3D(m_X_200);
  m_l_X_020 = global3DToLocal3D(m_X_020);
  m_l_X_002 = global3DToLocal3D(m_X_002);
  m_l_X_011 = global3DToLocal3D(m_X_011);
  m_l_X_101 = global3DToLocal3D(m_X_101);
  m_l_X_110 = global3DToLocal3D(m_X_110);
  
  m_coeff_x2 = m_l_X_020 - 2*m_l_X_110;
  m_coeff_y2 = m_l_X_002 - 2*m_l_X_101;
  m_coeff_xy = -2*m_l_X_110 + 2*m_l_X_011 - 2*m_l_X_101;
  m_coeff_x = 2*m_l_X_110;
  m_coeff_y = 2*m_l_X_101;
  
/*  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  qDebug()<<"m_X_200"<<m_X_200;
  qDebug()<<"m_X_020"<<m_X_020;
  qDebug()<<"m_X_002"<<m_X_002;
  qDebug()<<"m_X_011"<<m_X_011;
  qDebug()<<"m_X_101"<<m_X_101;
  qDebug()<<"m_X_110"<<m_X_110;
  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  qDebug()<<"m_l_X_200"<<m_l_X_200;
  qDebug()<<"m_l_X_020"<<m_l_X_020;
  qDebug()<<"m_l_X_002"<<m_l_X_002;
  qDebug()<<"m_l_X_011"<<m_l_X_011;
  qDebug()<<"m_l_X_101"<<m_l_X_101;
  qDebug()<<"m_l_X_110"<<m_l_X_110;
  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  qDebug()<<"m_coeff_x2"<<m_coeff_x2;
  qDebug()<<"m_coeff_y2"<<m_coeff_y2;
  qDebug()<<"m_coeff_xy"<<m_coeff_xy;
  qDebug()<<"m_coeff_x"<<m_coeff_x;
  qDebug()<<"m_coeff_y"<<m_coeff_y;
  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";*/
}

vec2_t BezierTriangle::fixedPointFunction(vec2_t t_inputPoint, double x, double y)
{
/*  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  qDebug()<<"m_coeff_x2"<<m_coeff_x2;
  qDebug()<<"m_coeff_y2"<<m_coeff_y2;
  qDebug()<<"m_coeff_xy"<<m_coeff_xy;
  qDebug()<<"m_coeff_x"<<m_coeff_x;
  qDebug()<<"m_coeff_y"<<m_coeff_y;
  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";*/
  vec2_t F;
  F[0] = pow(x,2)*m_coeff_x2[0] + pow(y,2)*m_coeff_y2[0] + x*y*m_coeff_xy[0] + x*m_coeff_x[0] + y*m_coeff_y[0] - t_inputPoint[0];
  F[1] = pow(x,2)*m_coeff_x2[1] + pow(y,2)*m_coeff_y2[1] + x*y*m_coeff_xy[1] + x*m_coeff_x[1] + y*m_coeff_y[1] - t_inputPoint[1];
  return F;
}

vec2_t BezierTriangle::fixedPointFunction(vec2_t t_inputPoint, vec2_t A)
{
  return fixedPointFunction(t_inputPoint, A[0], A[1]);
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

// mat3_t BezierTriangle::jacobiMatrix_no_projection(double x, double y)
// {
//   mat3_t J;
//   for(int i=0;i<3;i++) {
//     J[i][0] = 2*x*m_coeff_x2[i] + y*m_coeff_xy[i] + m_coeff_x[i];
//     J[i][1] = 2*y*m_coeff_y2[i] + x*m_coeff_xy[i] + m_coeff_y[i];
//   }
//   return J;
// }

mat2_t BezierTriangle::jacobiMatrix_numeric(vec2_t t_inputPoint, double x, double y, double dx, double dy)
{
  mat2_t J;
  vec2_t df = fixedPointFunction(t_inputPoint,x+dx,y+dy)-fixedPointFunction(t_inputPoint,x,y);
  if(dx<10e-9) {
    J[0][0] = 0;
    J[1][0] = 0;
  }
  else {
    J[0][0] = df[0]/dx;
    J[1][0] = df[1]/dx;
  }
  if(dy<10e-9) {
    J[0][1] = 0;
    J[1][1] = 0;
  }
  else {
    J[0][1] = df[0]/dy;
    J[1][1] = df[1]/dy;
  }
  return J;
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle3(vec3_t g_M)
{
  bool DEBUG = false;
  vec2_t t_M = global3DToLocal2D(g_M);
  if(DEBUG) qDebug()<<"t_M="<<t_M;
  vec2_t t_X = t_M;//vec2_t(1./3.,1./3.,1./3.);
  if(t_M[0]<0 || 1<t_M[0] || t_M[1]<0 || 1<t_M[1] || t_M[0]+t_M[1]>=1) {
    qDebug()<<"WARNING: Not on triangle! t_M="<<t_M;
    //TODO:
    //get closest point M' on triangle
    //get normal vector N at that point
    //project original point M onto plane (M',N)
  }
  
  if(DEBUG) qDebug()<<"t_X="<<t_X;
  vec2_t F = fixedPointFunction(t_M, t_X[0], t_X[1]);
  if(DEBUG) qDebug()<<"F.abs()="<<F.abs();
  int maxloops = 100;
  int Nloops=0;
  while(F.abs()>0.001 && Nloops < maxloops) {
    if(DEBUG) qDebug()<<"test passed with F.abs()="<<F.abs()<<" and "<<Nloops<<"<"<<maxloops;
    mat2_t J = jacobiMatrix(t_X[0], t_X[1]);
    if(J.det()==0) {
      qDebug()<<"WARNING: Matrix not invertible!";
    }
    if (fabs(J[0][0])+fabs(J[0][1])>=1) {
      if(DEBUG) qDebug()<<"WARNING: will not converge (case 1)";
    }
    if (fabs(J[1][0])+fabs(J[1][1])>=1) {
      if(DEBUG) qDebug()<<"WARNING: will not converge (case 2)";
    }
    
    mat2_t JI = J.inverse();
    vec2_t deltaX = -1*(JI*F);
    t_X = t_X + deltaX;
    if(DEBUG) qDebug()<<"t_X="<<t_X;
    F = fixedPointFunction(t_M, t_X[0], t_X[1]);
    if(DEBUG) qDebug()<<"F="<<F;
    if(DEBUG) qDebug()<<"F.abs()="<<F.abs();
    Nloops++;
  }
  if(Nloops >= maxloops) qDebug()<<"WARNING: Exited before converging! Nloops="<<Nloops;
  return QuadraticBezierTriangle(t_X);
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle4(vec3_t g_M)
{
  vec2_t t_M = global3DToLocal2D(g_M);
  qDebug()<<"t_M="<<t_M;

  vec2_t A,B,C,D,E,F;
  vec2_t K0,K1,K2;
  vec2_t delta;
  vec2_t x1,x2;
  vec2_t y1,y2;
  double x,y;
  
  A[0] = m_coeff_x2[0];
  B[0] = m_coeff_y2[0];
  C[0] = m_coeff_xy[0];
  D[0] = m_coeff_x[0];
  E[0] = m_coeff_y[0];
  F[0] = -t_M[0];
  
  A[1] = m_coeff_x2[1];
  B[1] = m_coeff_y2[1];
  C[1] = m_coeff_xy[1];
  D[1] = m_coeff_x[1];
  E[1] = m_coeff_y[1];
  F[1] = -t_M[1];
  
  vec2_t direction = fixedPointFunction(t_M, t_M[0], t_M[1]);
  if(direction[0]!=0) { // y=ax+b
    double a = direction[1]/direction[0];
    double b = -a*t_M[0];
    K0 = pow(b,2)*B + b*E + F;
    K1 = 2*a*b*B + b*C + D + a*E;
    K2 = A + pow(a,2)*B + a*C;
    // zero X component
    if(K2[0]!=0) {
      delta[0] = pow(K1[0],2)-4*K2[0]*K0[0];
      if(delta[0]>=0) {
        x1[0]=(-K1[0]+sqrt(delta[0]))/(2*K2[0]);
        x2[0]=(-K1[0]-sqrt(delta[0]))/(2*K2[0]);
      }
      else {
        qDebug()<<"FATAL ERROR: delta[0]<0";
        return vec3_t(0,0,0);
      }
    }
    else {
      x1[0]=-K0[0]/K1[0];
      x2[0]=x1[0];
    }
    // zero Y component
    if(K2[1]!=0) {
      delta[1] = pow(K1[1],2)-4*K2[1]*K0[1];
      if(delta[1]>=0) {
        x1[1]=(-K1[1]+sqrt(delta[1]))/(2*K2[1]);
        x2[1]=(-K1[1]-sqrt(delta[1]))/(2*K2[1]);
      }
      else {
        qDebug()<<"FATAL ERROR: delta[1]<0";
        return vec3_t(0,0,0);
      }
    }
    else {
      x1[1]=-K0[1]/K1[1];
      x2[1]=x1[1];
    }
    // calculate corresponding y values
    y1 = a*x1 + b*vec2_t(1,1);
    y2 = a*x2 + b*vec2_t(1,1);
  }
  else { // x=cte=ay+b=b=x0
    double a = 0;
    double b = t_M[0];
    K0 = pow(b,2)*A + b*D + F;
    K1 = b*C + E;
    K2 = B;
    // zero X component
    if(K2[0]!=0) {
      delta[0] = pow(K1[0],2)-4*K2[0]*K0[0];
      if(delta[0]>=0) {
        y1[0]=(-K1[0]+sqrt(delta[0]))/(2*K2[0]);
        y2[0]=(-K1[0]-sqrt(delta[0]))/(2*K2[0]);
      }
      else {
        qDebug()<<"FATAL ERROR: delta[0]<0";
        return vec3_t(0,0,0);
      }
    }
    else {
      y1[0]=-K0[0]/K1[0];
      y2[0]=y1[0];
    }
    // zero Y component
    if(K2[1]!=0) {
      delta[1] = pow(K1[1],2)-4*K2[1]*K0[1];
      if(delta[1]>=0) {
        y1[1]=(-K1[1]+sqrt(delta[1]))/(2*K2[1]);
        y2[1]=(-K1[1]-sqrt(delta[1]))/(2*K2[1]);
      }
      else {
        qDebug()<<"FATAL ERROR: delta[1]<0";
        return vec3_t(0,0,0);
      }
    }
    else {
      y1[1]=-K0[1]/K1[1];
      y2[1]=y1[1];
    }
    // calculate corresponding y values
    x1 = vec2_t(b,b);
    x2 = vec2_t(b,b);
  }
  
  qDebug()<<"x1="<<x1;
  qDebug()<<"x2="<<x2;
  qDebug()<<"y1="<<y1;
  qDebug()<<"y2="<<y2;
  
  vec2_t diff_1_00 = fixedPointFunction(t_M, x1[0], y1[0]);
  vec2_t diff_1_10 = fixedPointFunction(t_M, x1[1], y1[0]);
  vec2_t diff_1_01 = fixedPointFunction(t_M, x1[0], y1[1]);
  vec2_t diff_1_11 = fixedPointFunction(t_M, x1[1], y1[1]);
  
  vec2_t diff_2_00 = fixedPointFunction(t_M, x2[0], y2[0]);
  vec2_t diff_2_10 = fixedPointFunction(t_M, x2[1], y2[0]);
  vec2_t diff_2_01 = fixedPointFunction(t_M, x2[0], y2[1]);
  vec2_t diff_2_11 = fixedPointFunction(t_M, x2[1], y2[1]);
  
  qDebug()<<"diff_1_00.abs()="<<diff_1_00.abs();
  qDebug()<<"diff_1_10.abs()="<<diff_1_10.abs();
  qDebug()<<"diff_1_01.abs()="<<diff_1_01.abs();
  qDebug()<<"diff_1_11.abs()="<<diff_1_11.abs();
  
  qDebug()<<"diff_2_00.abs()="<<diff_2_00.abs();
  qDebug()<<"diff_2_10.abs()="<<diff_2_10.abs();
  qDebug()<<"diff_2_01.abs()="<<diff_2_01.abs();
  qDebug()<<"diff_2_11.abs()="<<diff_2_11.abs();
  
  vec2_t t_X(x,y);
  qDebug()<<"t_X="<<t_X;
  vec2_t diff = fixedPointFunction(t_M, t_X[0], t_X[1]);
  qDebug()<<"diff="<<diff;
  qDebug()<<"diff.abs()="<<diff.abs();
  return QuadraticBezierTriangle(t_X);
}

vec2_t secondDegreeSolver(double a, double b, double c) {
  double x1,x2;
  
  if(a==0) {
    if(b==0) {
      EG_BUG;
    }
    else {
      x1 = -c/b;
      x2 = x1;
    }
  }
  else {
    double delta = pow(b,2)-4*a*c;
    if(delta<0) {
      EG_BUG;
    }
    else {
      x1 = (-b+sqrt(delta))/(2*a);
      x2 = (-b-sqrt(delta))/(2*a);
    }
  }
  return vec2_t(x1,x2);
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle5(vec3_t g_M)
{
  vec2_t t_M = global3DToLocal2D(g_M);
  qDebug()<<"t_M="<<t_M;
  double x0 = t_M[0];
  double y0 = t_M[1];
  if(y0!=1) {
    vec2_t I(-x0/(y0-1),0);
//     vec2_t sol1 = secondDegreeSolver(A[0]+B[0]-2D[0],-2B[0]+2D[0],B[0]);
  }
  else {
    EG_BUG;
  }
  return vec3_t(0,0,0);
}

double BezierTriangle::z_func(vec2_t t_M)
{
  return z_func(t_M[0],t_M[1]);
}

double BezierTriangle::z_func(double x, double y)
{
  bool DEBUG = false;
  vec2_t t_M = vec2_t(x,y);
  if(DEBUG) qDebug()<<"t_M="<<t_M;
  vec2_t t_X = t_M;
  
  if(DEBUG) qDebug()<<"t_X="<<t_X;
  vec2_t F = fixedPointFunction(t_M, t_X[0], t_X[1]);
  if(DEBUG) qDebug()<<"F.abs()="<<F.abs();
  int maxloops = 100;
  int Nloops=0;
  while(F.abs()>0.001 && Nloops < maxloops) {
    if(DEBUG) qDebug()<<"test passed with F.abs()="<<F.abs()<<" and "<<Nloops<<"<"<<maxloops;
    mat2_t J = jacobiMatrix(t_X[0], t_X[1]);
    if(J.det()==0) {
      qDebug()<<"WARNING: Matrix not invertible!";
    }
    if (fabs(J[0][0])+fabs(J[0][1])>=1) {
      if(DEBUG) qDebug()<<"WARNING: will not converge (case 1)";
    }
    if (fabs(J[1][0])+fabs(J[1][1])>=1) {
      if(DEBUG) qDebug()<<"WARNING: will not converge (case 2)";
    }
    
    mat2_t JI = J.inverse();
    vec2_t deltaX = -1*(JI*F);
    t_X = t_X + deltaX;
    if(DEBUG) qDebug()<<"t_X="<<t_X;
    F = fixedPointFunction(t_M, t_X[0], t_X[1]);
    if(DEBUG) qDebug()<<"F="<<F;
    if(DEBUG) qDebug()<<"F.abs()="<<F.abs();
    Nloops++;
  }
  if(Nloops >= maxloops) qDebug()<<"WARNING: Exited before converging! Nloops="<<Nloops;
  vec3_t g_B = QuadraticBezierTriangle(t_X);
  vec3_t l_B = global3DToLocal3D(g_B);
  return l_B[2];
}

bool BezierTriangle::isInsideTriangle(vec2_t t_M, double tol)
{
  if(t_M[0]<0-tol || 1+tol<t_M[0] || t_M[1]<0-tol || 1+tol<t_M[1] || t_M[0]+t_M[1]>1+tol) {
    return false;
  }
  else {
    return true;
  }
}

vec3_t BezierTriangle::surfaceNormal(vec2_t t_M, int output)
{
  vec3_t bary_coords = getBarycentricCoordinates(t_M[0],t_M[1]);
  double u = bary_coords[0];
  double v = bary_coords[1];
  double w = bary_coords[2];
  
  vec2_t dx,dy;
  double k = 0.1 * smallest_length;
  vec2_t ex(1,0);
  vec2_t ey(0,1);
  if(u>=v && u>w) {
    dx = k * ex;
    dy = k * ey;
  }
  else if(v>=u && v>w) {
    dx = k * (ey-ex);
    dy = k * (-1*ex);
  }
  else if(w>=u && w>v) {
    dx = k * (-1*ey);
    dy = k * (ex-ey);
  }
  else {
    EG_BUG;
  }
  
/*  qDebug()<<"##############";
  qDebug()<<"dx="<<dx;
  qDebug()<<"dy="<<dy;
  qDebug()<<"##############";*/
  
  vec2_t t_P0 = t_M;
  double z0 = z_func(t_P0);
  vec3_t l_P0(t_P0[0], t_P0[1], z0 );
  
  vec2_t t_Px1 = t_P0-dx;
  vec2_t t_Px2 = t_P0+dx;
  vec2_t t_Py1 = t_P0-dy;
  vec2_t t_Py2 = t_P0+dy;
  double zx1 = z_func(t_Px1);
  double zx2 = z_func(t_Px2);
  double zy1 = z_func(t_Py1);
  double zy2 = z_func(t_Py2);
  vec3_t l_Px1(t_Px1[0],t_Px1[1],zx1);
  vec3_t l_Px2(t_Px2[0],t_Px2[1],zx2);
  vec3_t l_Py1(t_Py1[0],t_Py1[1],zy1);
  vec3_t l_Py2(t_Py2[0],t_Py2[1],zy2);
  
  vec3_t l_u1;
  vec3_t l_u2;
  
  if(!isInsideTriangle(t_P0)) {
    qWarning()<<"t_P0="<<t_P0;
    qWarning()<<"t_Px1="<<t_Px1;
    qWarning()<<"t_Px2="<<t_Px2;
    qWarning()<<"t_Py1="<<t_Py1;
    qWarning()<<"t_Py2="<<t_Py2;
    return vec3_t(0,0,0);
    EG_BUG;
  }
  
  if(isInsideTriangle(t_Px1) && isInsideTriangle(t_Px2)) {
    l_u1 = l_Px2-l_Px1;
  }
  else if(!isInsideTriangle(t_Px1) && isInsideTriangle(t_Px2)) {
    l_u1 = l_Px2-l_P0;
  }
  else if(isInsideTriangle(t_Px1) && !isInsideTriangle(t_Px2)) {
    l_u1 = l_P0-l_Px1;
  }
  else {
    qWarning()<<"t_P0="<<t_P0;
    qWarning()<<"t_Px1="<<t_Px1;
    qWarning()<<"t_Px2="<<t_Px2;
    qWarning()<<"t_Py1="<<t_Py1;
    qWarning()<<"t_Py2="<<t_Py2;
    return vec3_t(0,0,0);
    EG_BUG;
  }
  
  
  if(isInsideTriangle(t_Py1) && isInsideTriangle(t_Py2)) {
    l_u2 = l_Py2-l_Py1;
  }
  else if(!isInsideTriangle(t_Py1) && isInsideTriangle(t_Py2)) {
    l_u2 = l_Py2-l_P0;
  }
  else if(isInsideTriangle(t_Py1) && !isInsideTriangle(t_Py2)) {
    l_u2 = l_P0-l_Py1;
  }
  else {
    qWarning()<<"t_P0="<<t_P0;
    qWarning()<<"t_Px1="<<t_Px1;
    qWarning()<<"t_Px2="<<t_Px2;
    qWarning()<<"t_Py1="<<t_Py1;
    qWarning()<<"t_Py2="<<t_Py2;
    return vec3_t(0,0,0);
    EG_BUG;
  }
  
  vec3_t l_N = l_u1.cross(l_u2);
  vec3_t g_N = G*l_N;
  g_N.normalise();
  if(output==0) {
    return g_N;
  }
  else if(output==1) {
    return (G*l_u1).normalise();
  }
  else {
    return (G*l_u2).normalise();
  }
}
