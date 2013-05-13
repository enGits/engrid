// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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

void CadInterface::setForegroundGrid(vtkUnstructuredGrid *grid)
{
  m_FGrid = grid;
  m_FPart.trackGrid(m_FGrid);
}



vec3_t CadInterface::snap(vec3_t x, bool correct_curvature)

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

  if (correct_curvature) {
    x_best = correctCurvature(x_best);
  }

  m_LastNormal = n_best;
  m_LastRadius = r_best;
  return x_best;
}

vec3_t CadInterface::snapNode(vtkIdType id_node, bool correct_curvature)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  return snap(x, correct_curvature);
}

vec3_t CadInterface::snapNode(vtkIdType, vec3_t x, bool correct_curvature)
{
  return snap(x, correct_curvature);
}



vec3_t CadInterface::project(vec3_t x, vec3_t v, bool strict_direction, bool correct_curvature)
{
  vec3_t x_proj = x;
  m_LastNormal = v;
  m_LastRadius = 1e10;

  vec3_t x_hit1, n_hit1, x_hit2, n_hit2;
  double r_hit1, r_hit2;
  CadInterface::HitType hit_type1, hit_type2;

  hit_type1 = shootRay(x, v, x_hit1, n_hit1, r_hit1);
  if (hit_type1 == CadInterface::Miss && !strict_direction) {
    v *= -1;
    hit_type1 = shootRay(x, v, x_hit1, n_hit1, r_hit1);
  }
  if (hit_type1 == CadInterface::Miss) {
    m_Failed = true;
    return x;
  }
  m_Failed = false;
  v *= -1;
  x_proj = x_hit1;
  m_LastNormal = n_hit1;
  m_LastRadius = r_hit1;
  if (!strict_direction) {
    hit_type2 = shootRay(x_hit1, v, x_hit2, n_hit2, r_hit2);
    if (hit_type2 != CadInterface::Miss) {
      if ((x - x_hit2).abs() < (x - x_hit1).abs()) {
        x_proj = x_hit2;
        m_LastNormal = n_hit2;
        m_LastRadius = r_hit2;
      }
    }
  }

  if (correct_curvature) {
    x_proj = correctCurvature(x_proj);
  }

  return x_proj;
}

vec3_t CadInterface::projectNode(vtkIdType, vec3_t x, vec3_t v, bool strict_direction, bool correct_curvature)
{
  return project(x, v, strict_direction, correct_curvature);
}

vec3_t CadInterface::projectNode(vtkIdType id_node, vec3_t x, bool strict_direction, bool correct_curvature)
{
  vec3_t v = m_FPart.globalNormal(id_node);
  if (!checkVector(v)) {
    cout << "vector defect (id_node=" << id_node << ")" << endl;
    return x;
  }
  return projectNode(id_node, x, v, strict_direction, correct_curvature);
}

vec3_t CadInterface::projectNode(vtkIdType id_node, bool strict_direction, bool correct_curvature)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  return projectNode(id_node, x, strict_direction, correct_curvature);
}

double CadInterface::getRadius(vtkIdType id_node)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  snapNode(id_node, x, false);
  return m_LastRadius;
}

