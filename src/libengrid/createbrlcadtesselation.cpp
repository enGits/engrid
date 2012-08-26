//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "createbrlcadtesselation.h"

CreateBrlCadTesselation::CreateBrlCadTesselation(QString file_name, QString object_name)
{
  setupBrlCad(file_name, object_name);
}

bool CreateBrlCadTesselation::shootRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out)
{
  v.normalise();
  double r_hit;
  vec3_t x_hit, n_hit;
  BrlCadInterface::HitType hit_type = BrlCadInterface::shootRay(x, v, x_hit, n_hit, r_hit);
  if (hit_type == BrlCadInterface::Miss || hit_type == BrlCadInterface::HitOut) {
    return false;
  }
  x_in = x_hit;
  n_in = n_hit;
  x = x_in;
  do {
    x += 1e-10*v;
    hit_type = BrlCadInterface::shootRay(x, v, x_hit, n_hit, r_hit);
    if (hit_type == BrlCadInterface::HitOut) {
      x_out = x_hit;
      n_out = n_hit;
    }
    x = x_hit;
  } while (hit_type == BrlCadInterface::HitIn);
  return true;
}
