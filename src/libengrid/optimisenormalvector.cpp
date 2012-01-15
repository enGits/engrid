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

OptimiseNormalVector::OptimiseNormalVector(vec3_t n)
{
  m_InitialNormal = n;
}

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
  double err = 0;
  err += sqr(n.abs() - 1);
  n.normalise();
  foreach (vec3_t nc, m_Constraints) {
    err += sqr(n*nc);
  }
  if (m_Faces.size() > 0) {
    double h = 0;
    foreach (vec3_t nf, m_Faces) {
      h += n*nf;
    }
    h /= m_Faces.size();
    foreach (vec3_t nf, m_Faces) {
      err += sqr(h - n*nf);
    }
  }
  return err;
}

vec3_t OptimiseNormalVector::operator()()
{
  vec3_t n_opt = optimise(m_InitialNormal);
  if (!checkVector(n_opt)) {
    n_opt = m_InitialNormal;
  }
  return n_opt;
}
