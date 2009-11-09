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
}

Triangle::Triangle(vec3_t a_a, vec3_t a_b, vec3_t a_c)
{
  a = a_a;
  b = a_b;
  c = a_c;

  g1 = b-a;
  g2 = c-a;
  g3 = g1.cross(g2);
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
