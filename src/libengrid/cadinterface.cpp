//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                     +
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
#include "cadinterface.h"

#include <QTextStream>
#include <QFile>

CadInterface::CadInterface()
{
  QFile file(":/resources/misc/raysphere.dat");
  file.open(QIODevice::ReadOnly);
  QTextStream s(&file);
  int N;
  s >> N;
  m_RayPoints.resize(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < 3; ++j) {
      s >> m_RayPoints[i][j];
    }
    m_RayPoints[i].normalise();
  }
  m_LastRadius = 1e99;
  m_LastNormal = vec3_t(0,0,0);
}

vec3_t CadInterface::snap(vec3_t x)
{
  double L_best = 1e99;
  double r_best;
  vec3_t x_best;
  vec3_t n_best;
  bool no_hit = true;
  foreach (vec3_t x_ray, m_RayPoints) {
    double r_hit;
    vec3_t x_hit;
    vec3_t n_hit;
    vec3_t v = x_ray;

    // standard ray
    if (shootRay(x, v, x_hit, n_hit, r_hit) != Miss) {
      double L = (x_hit - x).abs();
      if (L < L_best) {
        x_best = x_hit;
        n_best = n_hit;
        r_best = r_hit;
        L_best = L;
        no_hit = false;
      }
      v = n_hit;

      // shoot a "far-side return ray" ...
      if (shootRay(x_hit, -1*x_ray, x_hit, n_hit, r_hit)) {
        double L = (x_hit - x).abs();
        if (L < L_best) {
          x_best = x_hit;
          n_best = n_hit;
          r_best = r_hit;
          L_best = L;
          no_hit = false;
        }
      }

      // re-shoot standard ray with surface normal
      if (shootRay(x, v, x_hit, n_hit, r_hit) != Miss) {
        double L = (x_hit - x).abs();
        if (L < L_best) {
          x_best = x_hit;
          n_best = n_hit;
          r_best = r_hit;
          L_best = L;
          no_hit = false;
        }
      }

      // re-shoot standard ray with surface normal in opposite direction
      v *= -1;
      if (shootRay(x, v, x_hit, n_hit, r_hit) != Miss) {
        double L = (x_hit - x).abs();
        if (L < L_best) {
          x_best = x_hit;
          n_best = n_hit;
          r_best = r_hit;
          L_best = L;
          no_hit = false;
        }
      }

    }
  }
  if (no_hit) {
    EG_BUG;
  }
  m_LastNormal = n_best;
  m_LastRadius = r_best;
  return x_best;
}
