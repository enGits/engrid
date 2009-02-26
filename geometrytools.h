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
#ifndef geometrytools_H
#define geometrytools_H

#include "math/mathvector.h"
#include "math/smallsquarematrix.h"

#include <vtkUnstructuredGrid.h>

namespace GeometryTools
{

/** Converts radians to degrees */
double rad2deg( double rad );

/** Converts degrees to radians */
double deg2rad( double deg );

void rotate(vec3_t g1, vec3_t g2, vec3_t g3, vec3_t &b, double theta);

/** Rotates vector v around vector axis by an angle theta */
vec3_t rotate(vec3_t v, vec3_t axis, double theta);

vec3_t orthogonalVector(vec3_t v);

double intersection(vec3_t x_straight, vec3_t v_straight, 
                    vec3_t x_plane, vec3_t n_plane);

double intersection(vec3_t x_straight, vec3_t v_straight, 
                    vec3_t x_plane, vec3_t u_plane, vec3_t v_plane);

bool intersection (double &k1, double &k2, vec2_t r1, vec2_t u1, vec2_t r2, vec2_t u2);

void sliceTriangle(const vector<vec3_t> &Tin, vec3_t x, vec3_t n, vector<vector<vec3_t> > &Tout);

double tetraVol(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, bool neg = false);

double pyraVol(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, bool neg = false);

double prismVol(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, vec3_t x6, bool neg = false);

double hexaVol(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4, vec3_t x5, vec3_t x6, vec3_t x7, vec3_t x8,  bool neg = false);

double triArea(vec3_t x1, vec3_t x2, vec3_t x3);

double quadArea(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4);

vec3_t triNormal(vec3_t x1, vec3_t x2, vec3_t x3);

vec3_t quadNormal(vec3_t x1, vec3_t x2, vec3_t x3, vec3_t x4);

vec3_t triNormal(vtkUnstructuredGrid *grid, vtkIdType p1, vtkIdType p2, vtkIdType p3);

vec3_t quadNormal(vtkUnstructuredGrid *grid, vtkIdType p1, vtkIdType p2, vtkIdType p3, vtkIdType p4);

vec3_t cellNormal(vtkUnstructuredGrid *grid, vtkIdType i);

double cellVA(vtkUnstructuredGrid *grid, vtkIdType cellId, bool neg = false);
  
inline vec2_t turnRight(const vec2_t &v)
{
  vec2_t u;
  u[0] =  v[1];
  u[1] = -v[0];
  return u;
};

inline vec2_t turnLeft(const vec2_t &v)
{
  vec2_t u;
  u[0] = -v[1];
  u[1] =  v[0];
  return u;
};

/** return the angle w.r.t. another 3-vector */
double angle(const vec3_t & u, const vec3_t & v);

};

#endif





