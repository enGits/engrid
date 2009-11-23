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
#include <gsl/gsl_poly.h>

BezierTriangle::BezierTriangle() : Triangle(), EgVtkObject() {
}

BezierTriangle::BezierTriangle(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110)  : Triangle(X_200, X_020, X_002), EgVtkObject() {
  setControlPoints(X_200,  X_020,  X_002,  X_011,  X_101,  X_110);
}

void BezierTriangle::setControlPoints(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110) {
  m_X_200 = X_200;
  m_X_020 = X_020;
  m_X_002 = X_002;
  m_X_011 = X_011;
  m_X_101 = X_101;
  m_X_110 = X_110;
  setupFunctionVariables();
}

void BezierTriangle::getControlPoints(vec3_t& X_200, vec3_t& X_020, vec3_t& X_002, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110) {
  X_200 = m_X_200;
  X_020 = m_X_020;
  X_002 = m_X_002;
  X_011 = m_X_011;
  X_101 = m_X_101;
  X_110 = m_X_110;
}

void BezierTriangle::writeBezierSurface(QString filename, int N) {
  //qDebug()<<"writeBezierSurface called";
//   int N=10;
  int N_cells = (N - 1) * (N - 1);
  int N_points = (N * N + N) / 2;

  //qDebug()<<"N_cells="<<N_cells;
  //qDebug()<<"N_points="<<N_points;

  EG_VTKSP(vtkUnstructuredGrid, bezier);
  allocateGrid(bezier, N_cells, N_points);

  vtkIdType offset = 0;
  offset += addBezierSurface(this, bezier, offset, N);

//   BezierTriangle B(m_X_200, m_X_020, m_X_002, m_X_011-vec3_t(0,0,1), m_X_101-vec3_t(0,0,1), m_X_110-vec3_t(0,0,1));
//   offset += B.addBezierSurface(bezier, offset, N);

  //qDebug()<<"offset="<<offset;

  saveGrid(bezier, filename);
}

vec3_t BezierTriangle::quadraticBezierTriangle(double u, double v, double w) {
  double total = u + v + w;
  u = u / total;
  v = v / total;
  w = w / total;
  return pow(u, 2)*m_X_200 + pow(v, 2)*m_X_020 + pow(w, 2)*m_X_002 + 2*u*v*m_X_110 + 2*v*w*m_X_011 + 2*w*u*m_X_101;
}

vec3_t BezierTriangle::quadraticBezierTriangle(vec2_t M) {
  vec3_t bary_coords = getBarycentricCoordinates(M[0], M[1]);
  double u, v, w;
  u = bary_coords[0];
  v = bary_coords[1];
  w = bary_coords[2];
  return quadraticBezierTriangle(u, v, w);
}

vec3_t BezierTriangle::quadraticBezierTriangle_g(vec3_t g_M) {
  vec2_t t_M = global3DToLocal2D(g_M);
  return quadraticBezierTriangle(t_M);
}

void BezierTriangle::setupFunctionVariables() {
  m_l_X_200 = global3DToLocal3D(m_X_200);
  m_l_X_020 = global3DToLocal3D(m_X_020);
  m_l_X_002 = global3DToLocal3D(m_X_002);
  m_l_X_011 = global3DToLocal3D(m_X_011);
  m_l_X_101 = global3DToLocal3D(m_X_101);
  m_l_X_110 = global3DToLocal3D(m_X_110);

  m_coeff_x2 = m_l_X_020 - 2 * m_l_X_110;
  m_coeff_y2 = m_l_X_002 - 2 * m_l_X_101;
  m_coeff_xy = -2 * m_l_X_110 + 2 * m_l_X_011 - 2 * m_l_X_101;
  m_coeff_x = 2 * m_l_X_110;
  m_coeff_y = 2 * m_l_X_101;

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

vec2_t BezierTriangle::fixedPointFunction(vec2_t t_inputPoint, double x, double y) {
  /*  qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
    qDebug()<<"m_coeff_x2"<<m_coeff_x2;
    qDebug()<<"m_coeff_y2"<<m_coeff_y2;
    qDebug()<<"m_coeff_xy"<<m_coeff_xy;
    qDebug()<<"m_coeff_x"<<m_coeff_x;
    qDebug()<<"m_coeff_y"<<m_coeff_y;
    qDebug()<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@";*/
  vec2_t F;
  F[0] = pow(x, 2) * m_coeff_x2[0] + pow(y, 2) * m_coeff_y2[0] + x * y * m_coeff_xy[0] + x * m_coeff_x[0] + y * m_coeff_y[0] - t_inputPoint[0];
  F[1] = pow(x, 2) * m_coeff_x2[1] + pow(y, 2) * m_coeff_y2[1] + x * y * m_coeff_xy[1] + x * m_coeff_x[1] + y * m_coeff_y[1] - t_inputPoint[1];
  return F;
}

vec2_t BezierTriangle::fixedPointFunction(vec2_t t_inputPoint, vec2_t A) {
  return fixedPointFunction(t_inputPoint, A[0], A[1]);
}

mat2_t BezierTriangle::jacobiMatrix(double x, double y) {
  mat2_t J;
  J[0][0] = 2 * x * m_coeff_x2[0] + y * m_coeff_xy[0] + m_coeff_x[0];
  J[1][0] = 2 * x * m_coeff_x2[1] + y * m_coeff_xy[1] + m_coeff_x[1];
  J[0][1] = 2 * y * m_coeff_y2[0] + x * m_coeff_xy[0] + m_coeff_y[0];
  J[1][1] = 2 * y * m_coeff_y2[1] + x * m_coeff_xy[1] + m_coeff_y[1];
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

mat2_t BezierTriangle::jacobiMatrix_numeric(vec2_t t_inputPoint, double x, double y, double dx, double dy) {
  mat2_t J;
  vec2_t df = fixedPointFunction(t_inputPoint, x + dx, y + dy) - fixedPointFunction(t_inputPoint, x, y);
  if (dx < 10e-9) {
    J[0][0] = 0;
    J[1][0] = 0;
  } else {
    J[0][0] = df[0] / dx;
    J[1][0] = df[1] / dx;
  }
  if (dy < 10e-9) {
    J[0][1] = 0;
    J[1][1] = 0;
  } else {
    J[0][1] = df[0] / dy;
    J[1][1] = df[1] / dy;
  }
  return J;
}

vec3_t BezierTriangle::projectLocal2DOnQuadraticBezierTriangle(vec2_t t_M) {
  bool DEBUG = false;

  vec2_t t_X = t_M;//vec2_t(1./3.,1./3.,1./3.);
  if (DEBUG) qDebug() << "t_X=" << t_X;
  vec2_t F = fixedPointFunction(t_M, t_X[0], t_X[1]);
  if (DEBUG) qDebug() << "F.abs()=" << F.abs();
  int maxloops = 100;
  int Nloops = 0;
  while (F.abs() > 0.001 && Nloops < maxloops) {
    if (DEBUG) qDebug() << "test passed with F.abs()=" << F.abs() << " and " << Nloops << "<" << maxloops;
    mat2_t J = jacobiMatrix(t_X[0], t_X[1]);
    if (J.det() == 0) {
      qDebug() << "WARNING: Matrix not invertible!";
    }
    if (fabs(J[0][0]) + fabs(J[0][1]) >= 1) {
      if (DEBUG) qDebug() << "WARNING: will not converge (case 1)";
    }
    if (fabs(J[1][0]) + fabs(J[1][1]) >= 1) {
      if (DEBUG) qDebug() << "WARNING: will not converge (case 2)";
    }

    mat2_t JI = J.inverse();
    vec2_t deltaX = -1 * (JI * F);
    t_X = t_X + deltaX;
    if (DEBUG) qDebug() << "t_X=" << t_X;
    F = fixedPointFunction(t_M, t_X[0], t_X[1]);
    if (DEBUG) qDebug() << "F=" << F;
    if (DEBUG) qDebug() << "F.abs()=" << F.abs();
    Nloops++;
  }
//   if (Nloops >= maxloops) qDebug() << "WARNING: Exited before converging! Nloops=" << Nloops;

  return quadraticBezierTriangle(t_X);
}

vec3_t BezierTriangle::projectOnQuadraticBezierTriangle(vec3_t g_M, int output) {
  vec2_t t_M = global3DToLocal2D(g_M);
//   qDebug()<<"t_M="<<t_M;

  if (!insideBezierSurface(g_M)) {
//     return vec3_t(0,0,0);
//     qDebug() << "WARNING: Not on bezier triangle! t_M=" << t_M;
    //get closest point M' on triangle
    int side = -1;
    double Lmin = 0;
    vec3_t g_Mp = closestPointOnBezierCurves(g_M, side, Lmin);
    
    if (!insideBezierSurface(g_Mp)) {
      setDebugLevel(1);
      insideBezierSurface(g_Mp);
      EG_BUG;
    }
    int zone = -1;
    vec3_t xi(0, 0, 0);
    vec3_t ri(0, 0, 0);
    double d = 0;
    projectOnTriangle(g_M, xi, ri, d, zone, true);
    if(zone==3 || zone==4 || zone==5) {
      side = zone;
    }
    
    /*    vec2_t t_Mp(ri[0], ri[1]);
    qDebug() << "t_Mp=" << t_Mp;
    vec3_t g_Mp = local2DToGlobal3D(t_Mp);*/
//     qDebug() << "g_Mp=" << g_Mp;
    vec2_t t_Mp = global3DToLocal2D(g_Mp);
    if(!checkVector(g_Mp)) EG_BUG;
    if(!checkVector(t_Mp)) EG_BUG;
    
//     if (output == 0) return g_Mp;
//     else return surfaceNormal(t_Mp, 0);
    
//     vec3_t g_Mp_proj = projectLocal2DOnQuadraticBezierTriangle(t_Mp);
//     qDebug() << "g_Mp_proj=" << g_Mp_proj;

//     qDebug()<<"side="<<side;
    if (m_has_neighbour[side]) {
      // no extrapolation, restrict
      if (output == 0) return g_Mp;
      else return surfaceNormal(t_Mp, 0);
    } else {
      // extrapolate

      //get normal vector N at that point
      vec3_t g_N = surfaceNormal(t_Mp, 0);
//       qDebug() << "g_N=" << g_N;

      //project original point M onto plane (M',N)
      if(!checkVector(g_Mp)) EG_BUG;
      if(!checkVector(g_N)) EG_BUG;
      double k = intersection(g_M, m_g3, g_Mp, g_N);
//     vec3_t g_P = projectPointOnPlane(g_M, g_Mp_proj, g_N);
      vec3_t g_P = g_M + k * m_g3;
      
      if(isnan(k) || isinf(k)) {
        qWarning()<<"g_M="<<g_M;
        qWarning()<<"m_g3="<<m_g3;
        qWarning()<<"g_Mp="<<g_Mp;
        qWarning()<<"g_N="<<g_N;
        EG_BUG;
      }
      if(!checkVector(g_M)) EG_BUG;
      if(!checkVector(m_g3)) EG_BUG;
      if(!checkVector(g_P)) EG_BUG;
      
      if (output == 0) return g_P;
      else return g_N;
    }

  } else {
    if (output == 0) return projectLocal2DOnQuadraticBezierTriangle(t_M);
    else return surfaceNormal(t_M, 0);
  }

}

vec2_t secondDegreeSolver(double a, double b, double c) {
  double x1, x2;

  if (a == 0) {
    if (b == 0) {
      EG_BUG;
    } else {
      x1 = -c / b;
      x2 = x1;
    }
  } else {
    double delta = pow(b, 2) - 4 * a * c;
    if (delta < 0) {
      EG_BUG;
    } else {
      x1 = (-b + sqrt(delta)) / (2 * a);
      x2 = (-b - sqrt(delta)) / (2 * a);
    }
  }
  return vec2_t(x1, x2);
}

double BezierTriangle::z_func(vec2_t t_M) {
  return z_func(t_M[0], t_M[1]);
}

double BezierTriangle::z_func(double x, double y) {
  bool DEBUG = false;
  vec2_t t_M = vec2_t(x, y);
  if (DEBUG) qDebug() << "t_M=" << t_M;
  vec3_t g_B = projectLocal2DOnQuadraticBezierTriangle(t_M);
  vec3_t l_B = global3DToLocal3D(g_B);
  return l_B[2];
}

vec3_t BezierTriangle::surfaceNormal(vec2_t t_M, int output) {
  vec3_t bary_coords = getBarycentricCoordinates(t_M[0], t_M[1]);
  double u = bary_coords[0];
  double v = bary_coords[1];
  double w = bary_coords[2];

  vec2_t dx, dy;
  vec2_t ex(1, 0);
  vec2_t ey(0, 1);
  if (u >= v && u > w) {
    dx = ex;
    dy = ey;
  } else if (v > u && v >= w) {
    dx = ey - ex;
    dy = -1 * ex;
  } else if (w >= u && w > v) {
    dx = -1 * ey;
    dy = ex - ey;
  } else {
    qWarning() << "bary_coords=" << bary_coords;
    EG_BUG;
  }

  /*  qDebug()<<"##############";
    qDebug()<<"dx="<<dx;
    qDebug()<<"dy="<<dy;
    qDebug()<<"##############";*/

  if (!isInsideTriangle(t_M)) {
    //TODO: special dx,dy for points outside basic triangle
    // get closest point on beziercurves
    int side = -1;
    double Lmin = 0;
    vec3_t g_M = local2DToGlobal3D(t_M);
    vec3_t g_Mp = closestPointOnBezierCurves(g_M, side, Lmin);
    vec2_t t_Mp = global3DToLocal2D(g_Mp);
    // get tangent at that point
    vec2_t t_tangent;
    insideBezierCurve(t_M, side, t_tangent);
    vec2_t t_normal = turnRight(t_tangent);// t_M-t_Mp;
//     vec2_t t_tangent = turnLeft(t_normal);
    checkVector(t_normal);
    checkVector(t_tangent);
    t_normal.normalise();
    t_tangent.normalise();
    checkVector(t_normal);
    checkVector(t_tangent);
    // build dx,dy according to that tangent + orig point + closest point
    dx = t_normal - t_tangent;
    dy = t_normal + t_tangent;
  }
  
  dx.normalise();
  dy.normalise();
  double k = 0.01;// * m_smallest_length;
  dx = k*dx;
  dy = k*dy;
  
  vec2_t t_P0 = t_M;
  double z0 = z_func(t_P0);
  vec3_t l_P0(t_P0[0], t_P0[1], z0);

  if (!insideBezierSurface(t_P0)) {
    qWarning() << "t_P0=" << t_P0;
  //     return vec3_t(0, 0, 0);
    EG_BUG;
  }
  
  vec3_t l_u1;
  vec3_t l_u2;
  
  int Nloops = 0;
  int maxloops = 100;
  bool dxdy_ok = false;
  while(!dxdy_ok && Nloops<maxloops) {
    
    dxdy_ok = true;
    
    vec2_t t_Px1 = t_P0 - dx;
    vec2_t t_Px2 = t_P0 + dx;
    vec2_t t_Py1 = t_P0 - dy;
    vec2_t t_Py2 = t_P0 + dy;
    double zx1 = z_func(t_Px1);
    double zx2 = z_func(t_Px2);
    double zy1 = z_func(t_Py1);
    double zy2 = z_func(t_Py2);
    vec3_t l_Px1(t_Px1[0], t_Px1[1], zx1);
    vec3_t l_Px2(t_Px2[0], t_Px2[1], zx2);
    vec3_t l_Py1(t_Py1[0], t_Py1[1], zy1);
    vec3_t l_Py2(t_Py2[0], t_Py2[1], zy2);
      
    if (insideBezierSurface(t_Px1) && insideBezierSurface(t_Px2)) {
      l_u1 = l_Px2 - l_Px1;
    } else if (!insideBezierSurface(t_Px1) && insideBezierSurface(t_Px2)) {
      l_u1 = l_Px2 - l_P0;
    } else if (insideBezierSurface(t_Px1) && !insideBezierSurface(t_Px2)) {
      l_u1 = l_P0 - l_Px1;
    } else {
/*      qWarning() << "t_P0=" << t_P0;
      qWarning() << "t_Px1=" << t_Px1;
      qWarning() << "t_Px2=" << t_Px2;
      qWarning() << "t_Py1=" << t_Py1;
      qWarning() << "t_Py2=" << t_Py2;*/
      dxdy_ok = false;
//       return vec3_t(0, 0, 0);
    }
  
  
    if (insideBezierSurface(t_Py1) && insideBezierSurface(t_Py2)) {
      l_u2 = l_Py2 - l_Py1;
    } else if (!insideBezierSurface(t_Py1) && insideBezierSurface(t_Py2)) {
      l_u2 = l_Py2 - l_P0;
    } else if (insideBezierSurface(t_Py1) && !insideBezierSurface(t_Py2)) {
      l_u2 = l_P0 - l_Py1;
    } else {
/*      qWarning() << "##############";
      qWarning() << "dx="<<dx;
      qWarning() << "dy="<<dy;
      qWarning() << "angle(dx,dy)=" << rad2deg(acos((dx*dy)/(dx.abs()*dy.abs())));
      qWarning() << "##############";
      qWarning() << "t_P0=" << t_P0 << " : " <<insideBezierSurface(t_P0);
      qWarning() << "t_Px1=" << t_Px1 << " : " <<insideBezierSurface(t_Px1);
      qWarning() << "t_Px2=" << t_Px2 << " : " <<insideBezierSurface(t_Px2);
      qWarning() << "t_Py1=" << t_Py1 << " : " <<insideBezierSurface(t_Py1);
      qWarning() << "t_Py2=" << t_Py2 << " : " <<insideBezierSurface(t_Py2);*/
      dxdy_ok = false;
//       return vec3_t(0, 0, 0);
    }
    
    if(!dxdy_ok) {
      dx = 0.5*dx;
      dy = 0.5*dy;
    }
    
    Nloops++;
  }
  if(Nloops>=maxloops) EG_BUG;
  
  vec3_t g_u1 = m_G * l_u1;
  g_u1.normalise();
  vec3_t g_u2 = m_G * l_u2;
  g_u2.normalise();
  vec3_t g_N = g_u1.cross(g_u2);
  g_N.normalise();
  if (output == 0) {
    return g_N;
  } else if (output == 1) {
    return g_u1;
  } else {
    return g_u2;
  }
}

bool BezierTriangle::insideBezierSurface(vec3_t g_M)
{
  vec2_t t_M = global3DToLocal2D(g_M);
//   qDebug()<<"g_M="<<g_M;
//   qDebug()<<"t_M="<<t_M;
  vec3_t xi(0, 0, 0);
  vec3_t ri(0, 0, 0);
  double d = 0;
  int side;
  projectOnTriangle(g_M, xi, ri, d, side, true);
  if(side==3 || side==4 || side==5) {
    if(DebugLevel>0) qWarning()<<"side==3 || side==4 || side==5";
    return false;
  }
  else {
    vec2_t t_tangent;
    if(insideBezierCurve(t_M,0,t_tangent) && insideBezierCurve(t_M,1,t_tangent) && insideBezierCurve(t_M,2,t_tangent)) {
      if(DebugLevel>0) qWarning()<<"insideBezierCurves";
      return true;
    }
    else {
      if(DebugLevel>0) qWarning()<<"else";
      return false;
    }
  }
}

bool BezierTriangle::insideBezierSurface(vec2_t t_M) {
  return insideBezierSurface(local2DToGlobal3D(t_M));
}

//TODO: merge projectOnBezierSide + insideBezierCurve ?
vec3_t BezierTriangle::projectOnBezierSide(vec3_t g_M, int side, double& Lmin, double& u)
{
  vec3_t a,b,c;
  if(side==0) { // w=0
    // B-M = a*u^2 + b*u + c
    a = (m_X_200 + m_X_020 - 2*m_X_110);
    b = (-2*m_X_020 + 2*m_X_110);
    c = m_X_020 - g_M;
  }
  else if(side==1) { // u=0
    a = (m_X_020 + m_X_002 - 2*m_X_011);
    b = (-2*m_X_002 + 2*m_X_011);
    c = m_X_002 - g_M;
  }
  else { // v=0
    a = (m_X_002 + m_X_200 - 2*m_X_101);
    b = (-2*m_X_200 + 2*m_X_101);
    c = m_X_200 - g_M;
  }
  
  // d((B-M).abs())/du = coeff3*u^3 + coeff2*u^2 + coeff1*u + coeff0
  double coeff3, coeff2, coeff1, coeff0;
  coeff3 = 4*a*a;
  coeff2 = 6*a*b;
  coeff1 = 4*a*c+2*b*b;
  coeff0 = 2*b*c;
  
  double x[3];
  int N;
  if(coeff3!=0) N = gsl_poly_solve_cubic(coeff2/coeff3, coeff1/coeff3, coeff0/coeff3, &(x[0]), &(x[1]), &(x[2]));
  else N = gsl_poly_solve_quadratic (coeff2, coeff1, coeff0, &(x[0]), &(x[1]));
  
  if(N==0) EG_BUG;
  
  double L[3];
  Lmin = 0;
  u = 0;
  bool first = true;
  
  for(int i=0;i<N;i++) {
    if(x[i]<0) x[i]=0;
    if(x[i]>1) x[i]=1;
    L[i] = (pow(x[i],2)*a + x[i]*b + c).abs();
    if(first) {
      Lmin = L[i];
      u = x[i];
      first = false;
    }
    else if(L[i]<Lmin) {
      Lmin = L[i];
      u = x[i];
    }
  }
  
  vec3_t g_B = g_M + pow(u,2)*a + u*b + c;
  return g_B;
}

bool BezierTriangle::insideBezierCurve(vec2_t t_M, int side, vec2_t& t_tangent, double tol)
{
//   qDebug()<<"t_M="<<t_M;
//   qDebug()<<"side="<<side;
  
  vec2_t t_X_200 = global3DToLocal2D(m_X_200);
  vec2_t t_X_020 = global3DToLocal2D(m_X_020);
  vec2_t t_X_002 = global3DToLocal2D(m_X_002);
  vec2_t t_X_011 = global3DToLocal2D(m_X_011);
  vec2_t t_X_101 = global3DToLocal2D(m_X_101);
  vec2_t t_X_110 = global3DToLocal2D(m_X_110);
  
  vec2_t a,b,c;
  if(side==0) { // w=0
    // B -M = a*u^2 + b*u + c
    a = (t_X_200 + t_X_020 - 2*t_X_110);
    b = (-2*t_X_020 + 2*t_X_110);
    c = t_X_020 - t_M;
  }
  else if(side==1) { // u=0
    a = (t_X_020 + t_X_002 - 2*t_X_011);
    b = (-2*t_X_002 + 2*t_X_011);
    c = t_X_002 - t_M;
  }
  else { // v=0
    a = (t_X_002 + t_X_200 - 2*t_X_101);
    b = (-2*t_X_200 + 2*t_X_101);
    c = t_X_200 - t_M;
  }
  
  // d((B-M).abs())/du = coeff3*u^3 + coeff2*u^2 + coeff1*u + coeff0
  double coeff3, coeff2, coeff1, coeff0;
  coeff3 = 4*a*a;
  coeff2 = 6*a*b;
  coeff1 = 4*a*c+2*b*b;
  coeff0 = 2*b*c;
  
  double x[3];
  int N;
  if(coeff3!=0) N = gsl_poly_solve_cubic(coeff2/coeff3, coeff1/coeff3, coeff0/coeff3, &(x[0]), &(x[1]), &(x[2]));
  else N = gsl_poly_solve_quadratic (coeff2, coeff1, coeff0, &(x[0]), &(x[1]));
  
  if(N==0) EG_BUG;
  
  double L[3];
  double Lmin = 0;
  double u = 0;
  bool first = true;
  
  for(int i=0;i<N;i++) {
    if(isnan(x[i]) || isinf(x[i])) {
      qWarning()<<"NAN OR INF";
      qWarning()<<"x[i]="<<x[i];
      qWarning()<<"coeff3="<<coeff3;
      qWarning()<<"coeff2="<<coeff2;
      qWarning()<<"coeff1="<<coeff1;
      qWarning()<<"coeff0="<<coeff0;
      EG_BUG;
    }
    if(x[i]<0) x[i]=0;
    if(x[i]>1) x[i]=1;
    L[i] = (pow(x[i],2)*a + x[i]*b + c).abs();
    if(first) {
      Lmin = L[i];
      u = x[i];
      first = false;
    }
    else if(L[i]<Lmin) {
      Lmin = L[i];
      u = x[i];
    }
  }
  
  vec2_t t_B = t_M + pow(u,2)*a + u*b + c;
  t_tangent = 2*u*a + b;
  
  vec3_t l_M = vec3_t(t_M[0],t_M[1],0);
  vec3_t l_B = vec3_t(t_B[0],t_B[1],0);
  vec3_t l_tangent = vec3_t(t_tangent[0],t_tangent[1],0);
  
  checkVector(l_B);
  checkVector(l_tangent);
  if(DebugLevel>0) {
    qWarning()<<"l_M="<<l_M;
    qWarning()<<"l_B="<<l_B;
    qWarning()<<"l_tangent="<<l_tangent;
    qWarning()<<"l_tangent.cross(l_M-l_B)="<<l_tangent.cross(l_M-l_B);
    qWarning()<<"(l_tangent.cross(l_M-l_B))[2]="<<(l_tangent.cross(l_M-l_B))[2];
  }
  return ( (l_tangent.cross(l_M-l_B))[2]<=0+tol );
}

vec3_t BezierTriangle::closestPointOnBezierCurves(vec3_t g_M, int& side, double& Lmin)
{
  Lmin = 0;
  side = -1;
  vec3_t g_Mp;
  bool first = true;
  for(int i_side=0; i_side<3; i_side++) {
    double L,u;
    vec3_t foo = projectOnBezierSide(g_M, i_side,L,u);
    if(first || L<Lmin) {
      Lmin = L;
      g_Mp = foo;
      side = i_side;
      first = false;
    }
  }
  return g_Mp;
}
