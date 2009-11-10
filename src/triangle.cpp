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
#include "triangle.h"
#include "geometrytools.h"
#include "engrid.h"

Triangle::Triangle()
{
  id_a=0;
  id_b=0;
  id_c=0;
  a=vec3_t(0,0,0);
  b=vec3_t(0,0,0);
  c=vec3_t(0,0,0);
  setupTriangle();
}

Triangle::Triangle(vec3_t a_a, vec3_t a_b, vec3_t a_c)
{
  a = a_a;
  b = a_b;
  c = a_c;
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_a, vtkIdType a_id_b, vtkIdType a_id_c)
{
  id_a = a_id_a;
  id_b = a_id_b;
  id_c = a_id_c;
  a_grid->GetPoints()->GetPoint(id_a, a.data());
  a_grid->GetPoints()->GetPoint(id_b, b.data());
  a_grid->GetPoints()->GetPoint(id_c, c.data());
  setupTriangle();
}

Triangle::Triangle(vtkUnstructuredGrid* a_grid, vtkIdType a_id_cell)
{
  vtkIdType Npts, *pts;
  a_grid->GetCellPoints(a_id_cell, Npts, pts);
  if (Npts == 3) {
    id_a = pts[0];
    id_b = pts[1];
    id_c = pts[2];
    a_grid->GetPoints()->GetPoint(id_a, a.data());
    a_grid->GetPoints()->GetPoint(id_b, b.data());
    a_grid->GetPoints()->GetPoint(id_c, c.data());
    setupTriangle();
  } else {
    EG_ERR_RETURN("only triangles allowed at the moment");
  }
}

void Triangle::setupTriangle()
{
  g1 = b - a;
  g2 = c - a;
  g3 = g1.cross(g2);
  
  A  = 0.5*g3.abs();
  g3.normalise();
  
  G.column(0, g1);
  G.column(1, g2);
  G.column(2, g3);
  GI = G.inverse();
  
  smallest_length = (b - a).abs();
  smallest_length = min(smallest_length, (c - b).abs());
  smallest_length = min(smallest_length, (a - c).abs());
}

vec3_t Triangle::local3DToGlobal3D(vec3_t l_M)
{
  return a+G*l_M;
}

vec3_t Triangle::global3DToLocal3D(vec3_t g_M)
{
  vec3_t tmp = g_M-a;
  return GI*tmp;
}

vec3_t Triangle::local2DToGlobal3D(vec2_t l_M)
{
  return local3DToGlobal3D(vec3_t(l_M[0],l_M[1],0));
}

vec2_t Triangle::global3DToLocal2D(vec3_t g_M)
{
  vec3_t l_M = global3DToLocal3D(g_M);
  return vec2_t(l_M[0],l_M[1]);
}

bool Triangle::projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d)
{
  xi = vec3_t(1e99,1e99,1e99);
  double scal = (xp - this->a)*this->g3;
  vec3_t x1, x2;
  if (scal > 0) {
    x1 = xp + this->g3;
    x2 = xp - scal*this->g3 - this->g3;
  } else {
    x1 = xp - this->g3;
    x2 = xp - scal*this->g3 + this->g3;
  }
  d = 1e99;
  bool intersects_face = GeometryTools::intersectEdgeAndTriangle(this->a, this->b, this->c, x1, x2, xi, ri);
  if (intersects_face) {
    vec3_t dx = xp - this->a;
    d = fabs(dx*this->g3);
  } else {
    double kab = GeometryTools::intersection(this->a, this->b - this->a, xp, this->b - this->a);
    double kac = GeometryTools::intersection(this->a, this->c - this->a, xp, this->c - this->a);
    double kbc = GeometryTools::intersection(this->b, this->c - this->b, xp, this->c - this->b);
    double dab = (this->a + kab*(this->b-this->a) - xp).abs();
    double dac = (this->a + kac*(this->c-this->a) - xp).abs();
    double dbc = (this->b + kbc*(this->c-this->b) - xp).abs();
    bool set = false;
    if ((kab >= 0) && (kab <= 1)) {
      if (dab < d) {
        xi = this->a + kab*(this->b-this->a);
        d = dab;
        set = true;
      }
    }
    if ((kac >= 0) && (kac <= 1)) {
      if (dac < d) {
        xi = this->a + kac*(this->c-this->a);
        d = dac;
        set = true;
      }
    }
    if ((kbc >= 0) && (kbc <= 1)) {
      if (dbc < d) {
        xi = this->b + kbc*(this->c-this->b);
        d = dbc;
        set = true;
      }
    }
    double da = (this->a - xp).abs();
    double db = (this->b - xp).abs();
    double dc = (this->c - xp).abs();
    if (da < d) {
      xi = this->a;
      d = da;
      set = true;
    }
    if (db < d) {
      xi = this->b;
      d = db;
    }
    if (dc < d) {
      xi = this->c;
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
