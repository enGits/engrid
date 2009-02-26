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
#include "geometrytools.h"
#include "containertricks.h"

#include <vtkCellType.h>
#include <cmath>

namespace GeometryTools {

double rad2deg( double rad )
{
  return rad/M_PI*180;
}

double deg2rad( double deg )
{
  return deg/180*M_PI;
}

void rotate(vec3_t g1, vec3_t g2,vec3_t g3, vec3_t &b, double theta)
{
  //cout << b << endl;
  mat3_t g_t;
  g_t[0] = g1;
  g_t[1] = g2;
  g_t[2] = g3;
  mat3_t g = g_t.transp();
  //cout << g << endl;
  vec3_t bb = g_t*b;
  //cout << bb << endl;
  mat3_t rot = mat3_t::identity();
  rot[0][0] = cos(theta);
  rot[0][1] = -sin(theta);
  rot[1][0] = sin(theta);
  rot[1][1] = cos(theta);
  //cout << rot << endl;
  vec3_t bbb = rot*bb;
  //cout << bbb << endl;
  b = g*bbb;
  //cout << b << endl;
};

vec3_t rotate(vec3_t v, vec3_t axis, double theta)
{
  axis.normalise();
  
  // transposed base of rotate system
  mat3_t g_t;

  // compute projection of v in axis direction
  vec3_t v_axis = (axis*v)*axis;

  // compute first orthogonal vector (first base vector)
  g_t[0] = v-v_axis;
  
  //In case of points on the rotation axis, do nothing
  if(g_t[0].abs()==0) return v;
  
  g_t[0].normalise();

  // second base vector is the normalised axis
  g_t[1] = axis;

  // compute second orthogonal vector (third base vector)
  g_t[2] = g_t[0].cross(g_t[1]);

  // base of rotate system
  mat3_t g = g_t.transp();

  // matrix for rotation around g_t[1];
  mat3_t rot = mat3_t::identity();
  rot[0][0] = cos(theta);
  rot[0][2] = sin(theta);
  rot[2][0] = -sin(theta);
  rot[2][2] = cos(theta);

  // transfer v to rotate system
  vec3_t v_r = g_t*v;

  // rotate the vector and transfer it back
  v_r = rot*v_r;
  v = g*v_r;

  return v;
};

vec3_t orthogonalVector(vec3_t v)
{
  v.normalise();
  vec3_t u = v;
  int i_min = 0;
  int i_max = 0;
  for (int i = 1; i < 3; ++i) {
    if (v[i] > v[i_max]) i_max = i;
    if (v[i] < v[i_min]) i_min = i;
  };
  double h = u[i_min];
  u[i_min] = u[i_max];
  u[i_max] = h;
  u -= (u*v)*v;
  u.normalise();
  return u;
};


double intersection(vec3_t x_straight, vec3_t v_straight, vec3_t x_plane, vec3_t u_plane, vec3_t v_plane)
{
  vec3_t n = u_plane.cross(v_plane);
  return intersection(x_straight,v_straight,x_plane,n);
};

double intersection(vec3_t x_straight, vec3_t v_straight, vec3_t x_plane, vec3_t n_plane)
{
  double k = (x_plane*n_plane - x_straight*n_plane)/(v_straight*n_plane);
  return k;
};

bool intersection (double &k1, double &k2, vec2_t r1, vec2_t u1, vec2_t r2, vec2_t u2)
{
  double ave_length = .5*(sqrt(u1[0]*u1[0]+u1[1]*u1[1]) + sqrt(u2[0]*u2[0]+u2[1]*u2[1]));
  double DET = (u1[0]*u2[1]-u1[1]*u2[0]);
  if (fabs(DET) > 1e-6*ave_length) {
    k1 = -(u2[0]*r2[1]-u2[0]*r1[1]-r2[0]*u2[1]+r1[0]*u2[1])/DET;
    k2 = -(-u1[1]*r2[0]+u1[0]*r2[1]-u1[0]*r1[1]+u1[1]*r1[0])/DET;
    return true;
  } else {
    return false;
  };
};

void sliceTriangle(const vector<vec3_t> &Tin, vec3_t x, vec3_t n, vector<vector<vec3_t> > &Tout)
{
  vec3_t a = Tin[0];
  vec3_t b = Tin[1];
  vec3_t c = Tin[2];
  double kab = intersection(a,b-a,x,n);
  double kbc = intersection(b,c-b,x,n);
  double kca = intersection(c,a-c,x,n);
  bool ab_cut = ((kab >= 0) && (kab <= 1));
  bool bc_cut = ((kbc >= 0) && (kbc <= 1));
  bool ca_cut = ((kca >= 0) && (kca <= 1));
  if (ab_cut && bc_cut && ca_cut) {
    //cerr << "invalid triangle (SliceTriangle) A" << endl;
    //exit(EXIT_FAILURE);
    if      ((kab <= kbc) && (kab <= kca)) ab_cut = false;
    else if ((kbc <= kab) && (kbc <= kca)) bc_cut = false;
    else                                   ca_cut = false;
  };
  if (ab_cut && bc_cut) {
    vec3_t ab = a + kab*(b-a);
    vec3_t bc = b + kbc*(c-b);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,ab,bc;
    clinit(Tout[1]) = ab,b,bc;
    clinit(Tout[2]) = bc,c,a;
  } else if (bc_cut && ca_cut) {
    vec3_t bc = b + kbc*(c-b);
    vec3_t ca = c + kca*(a-c);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,bc,ca;
    clinit(Tout[1]) = a,b,bc;
    clinit(Tout[2]) = bc,c,ca;
  } else if (ca_cut && ab_cut) {
    vec3_t ca = c + kca*(a-c);
    vec3_t ab = a + kab*(b-a);
    Tout.resize(3,vector<vec3_t>(3));
    clinit(Tout[0]) = a,ab,ca;
    clinit(Tout[1]) = ab,b,ca;
    clinit(Tout[2]) = b,c,ca;
  } else {
    Tout.resize(1,vector<vec3_t>(3));
    clinit(Tout[0]) = a,b,c;
  };
};


double tetraVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, bool neg)
{
  vec3_t v1 = x1-x0;
  vec3_t v2 = x2-x0;
  vec3_t v3 = x3-x0;
  double V = v1*(v2.cross(v3));
  V /= 6.0;
  if (!neg && (V < 0)) {
    vec3_t v4 = x2-x1;
    vec3_t v5 = x3-x1;
    vec3_t v6 = x3-x2;
    double Lmin = 1e99;
    double Lmax = 0;
    Lmin = min(Lmin,v1.abs());
    Lmin = min(Lmin,v2.abs());
    Lmin = min(Lmin,v3.abs());
    Lmin = min(Lmin,v4.abs());
    Lmin = min(Lmin,v5.abs());
    Lmin = min(Lmin,v6.abs());
    Lmax = max(Lmax,v1.abs());
    Lmax = max(Lmax,v2.abs());
    Lmax = max(Lmax,v3.abs());
    Lmax = max(Lmax,v4.abs());
    Lmax = max(Lmax,v5.abs());
    Lmax = max(Lmax,v6.abs());
    if (Lmin/Lmax < 1e-6) {
      V = 0;
      V = -1e99;
    } else {
      V = -1e99;
    };
    //V *= 1e6;
  };
  return V; //fabs(V);
};

double pyraVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, bool neg)
{
  double V1 = tetraVol(x0, x1, x3, x4, neg) + tetraVol(x1, x2, x3, x4, neg);
  double V2 = tetraVol(x0, x1, x2, x4, neg) + tetraVol(x2, x3, x0, x4, neg);
  return min(V1,V2);
  /*
  double V = 0;
  vec3_t m0 = .25*(x0+x1+x2+x3);
  V += tetraVol(x0, x1, m0, x4, neg);
  V += tetraVol(x1, x2, m0, x4, neg);
  V += tetraVol(x2, x3, m0, x4, neg);
  V += tetraVol(x3, x0, m0, x4, neg);
  return V;
  */
};

double prismVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, bool neg)
{
  double V = 0;
  vec3_t p  = 1.0/6.0*(x0+x1+x2+x3+x4+x5);
  V += tetraVol(x0, x2, x1, p, neg);
  V += tetraVol(x3, x4, x5, p, neg);
  V += pyraVol (x0, x1, x4, x3, p, neg);
  V += pyraVol (x1, x2, x5, x4, p, neg);
  V += pyraVol (x0, x3, x5, x2, p, neg);
  return V;
};

double hexaVol(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, vec3_t x6, vec3_t x7, bool neg)
{
  double V = 0;
  vec3_t p = 1.0/8.0*(x0+x1+x2+x3+x4+x5+x6+x7);
  V += pyraVol(x0, x1, x3, x2, p, neg);
  V += pyraVol(x0, x4, x5, x1, p, neg);
  V += pyraVol(x4, x6, x7, x5, p, neg);
  V += pyraVol(x2, x3, x7, x6, p, neg);
  V += pyraVol(x1, x5, x7, x3, p, neg);
  V += pyraVol(x0, x2, x6, x4, p, neg);
  return V;
};

double triArea(vec3_t x0, vec3_t x1, vec3_t x2)
{
  vec3_t a = x1-x0;
  vec3_t b = x2-x0;
  double A = 0.5*((a.cross(b)).abs());
  return A;
};

double quadArea(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3)
{
  double A = 0;
  vec3_t p = .25*(x0+x1+x2+x3);
  A += triArea(x0,x1,p);
  A += triArea(x1,x2,p);
  A += triArea(x2,x3,p);
  A += triArea(x3,x0,p);
  return A;
};

vec3_t triNormal(vec3_t x0, vec3_t x1, vec3_t x2)
{
  vec3_t a = x1-x0;
  vec3_t b = x2-x0;
  vec3_t n = 0.5*(a.cross(b));
  return n;
};

vec3_t quadNormal(vec3_t x0, vec3_t x1, vec3_t x2, vec3_t x3)
{
  vec3_t n;
  clinit(n) = 0,0,0;
  vec3_t p = .25*(x0+x1+x2+x3);
  n += triNormal(x0,x1,p);
  n += triNormal(x1,x2,p);
  n += triNormal(x2,x3,p);
  n += triNormal(x3,x0,p);
  return n;
};

vec3_t triNormal(vtkUnstructuredGrid *grid, vtkIdType p1, vtkIdType p2, vtkIdType p3)
{
  vec3_t x1, x2, x3;
  grid->GetPoint(p1,x1.data());
  grid->GetPoint(p2,x2.data());
  grid->GetPoint(p3,x3.data());
  return triNormal(x1,x2,x3);
};

vec3_t quadNormal(vtkUnstructuredGrid *grid, vtkIdType p1, vtkIdType p2, vtkIdType p3, vtkIdType p4)
{
  vec3_t x1, x2, x3, x4;
  grid->GetPoint(p1,x1.data());
  grid->GetPoint(p2,x2.data());
  grid->GetPoint(p3,x3.data());
  grid->GetPoint(p4,x4.data());
  return quadNormal(x1,x2,x3,x4);
};

vec3_t cellNormal(vtkUnstructuredGrid *grid, vtkIdType i)
{
  vtkIdType *pts;
  vtkIdType npts;
  vec3_t n;
  clinit(n) = 0,0,0;
  grid->GetCellPoints(i,npts,pts);
  if      (npts == 3) return triNormal(grid,pts[0],pts[1],pts[2]);
  else if (npts == 4) return quadNormal(grid,pts[0],pts[1],pts[2],pts[3]);
  return n;
};

double cellVA(vtkUnstructuredGrid *grid, vtkIdType cellId, bool neg)
{
  vtkIdType *pts;
  vtkIdType  Npts;
  vec3_t     p[8];
  grid->GetCellPoints(cellId, Npts, pts);
  for (int i_pts = 0; i_pts < Npts; ++i_pts) {
    grid->GetPoints()->GetPoint(pts[i_pts], p[i_pts].data());
  };
  vtkIdType cellType = grid->GetCellType(cellId);
  if      (cellType == VTK_TRIANGLE)   return triArea (p[0], p[1], p[2]);
  else if (cellType == VTK_QUAD)       return quadArea(p[0], p[1], p[2], p[3]);
  else if (cellType == VTK_TETRA)      return tetraVol(p[0], p[1], p[2], p[3], neg);
  else if (cellType == VTK_PYRAMID)    return pyraVol (p[0], p[1], p[2], p[3], p[4], neg);
  else if (cellType == VTK_WEDGE)      return prismVol(p[0], p[1], p[2], p[3], p[4], p[5], neg);
  else if (cellType == VTK_HEXAHEDRON) return hexaVol (p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], neg);
  return 0.0;
};

double angle(const vec3_t & u, const vec3_t & v)
{
  // return the angle w.r.t. another 3-vector
  double ptot2 = u.abs2()*v.abs2();
  if(ptot2 <= 0) {
      return 0.0;
  } else {
    double arg = (u*v)/sqrt(ptot2);
    if(arg >  1.0) arg =  1.0;
    if(arg < -1.0) arg = -1.0;
    return acos(arg);
  }
}

} // namespace




