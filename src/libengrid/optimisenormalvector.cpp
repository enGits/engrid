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
#include "optimisenormalvector.h"

void OptimiseNormalVector::addConstraint(vec3_t n)
{
  m_Constraints.append(n);
}

void OptimiseNormalVector::addFace(vec3_t n)
{
  m_Faces.append(n);
}

double OptimiseNormalVector::func(vec3_t n)
{
  double hf = 1;
  double hc = 0;
  vec3_t n0 = n;
  n0.normalise();
  foreach (vec3_t nc, m_Constraints) {
    nc.normalise();
    double h = nc*n0;
    hc = max(h, fabs(hc));
  }
  foreach (vec3_t nf, m_Faces) {
    nf.normalise();
    double h = nf*n0;
    hf = min(h, hf);
  }
  return sqr(hc) + sqr(1 - hf) + sqr(1 - n.abs());
}

vec3_t OptimiseNormalVector::operator()(vec3_t n)
{
  vec3_t n_opt = optimise(n);
  if (!checkVector(n_opt)) {
    n_opt = n;
  }
  n_opt.normalise();
  return n_opt;
}
