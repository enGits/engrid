// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
  double hf_min =  1;
  double hf_max = -1;
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
    hf_min = min(h, hf_min);
    hf_max = max(h, hf_max);
  }
  return sqr(hc) + sqr(1 - hf_min) + sqr(1 - hf_max);
}

vec3_t OptimiseNormalVector::optimise(vec3_t n)
{
  int count = 0;
  computeDerivatives(n);
  n.normalise();
  double scale = 0.1;
  while (count < 100 && scale > 1e-4) {
    double ag = grad_f.abs();
    double err1 = func(n);
    vec3_t dn = -1.0*grad_f;
    dn -= (n*dn)*n;
    if (grad_f.abs() > 1e-10) {
      dn.normalise();
    }
    double relax = min(scale, scale*grad_f.abs());
    dn *= relax;
    vec3_t n_old = n;
    n += dn;
    n.normalise();
    double err2 = func(n);
    if (err2 > err1) {
      scale *= 0.1;
      n = n_old;
    }
    ++count;
    computeDerivatives(n);
  }
  return n;
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
