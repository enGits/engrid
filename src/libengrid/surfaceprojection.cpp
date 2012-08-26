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
#include "surfaceprojection.h"

SurfaceProjection::SurfaceProjection(CadInterface *cad_interface)
{
  m_CadInterface = cad_interface;
}

void SurfaceProjection::setForegroundGrid(vtkUnstructuredGrid *grid)
{
  m_FGrid = grid;
  m_FPart.trackGrid(m_FGrid);
}

vec3_t SurfaceProjection::project(vec3_t x, vtkIdType id_node, bool, vec3_t v, bool strict_direction, bool allow_search)
{
  vec3_t n = v;
  if (n.abs() < 1e-3) {
    if (id_node == -1) {
      EG_BUG;
    }
    n = m_FPart.globalNormal(id_node);
  }
  if (!checkVector(x)) {
    EG_BUG;
  }
  if (!checkVector(n)) {
    cout << "vector defect (id_node=" << id_node << ")" << endl;
    return x;
    EG_BUG;
  }

  vec3_t x_proj = x;
  m_LastNormal = n;
  m_LastRadius = 1e10;

  vec3_t x_hit1, n_hit1, x_hit2, n_hit2;
  double r_hit1, r_hit2;
  CadInterface::HitType hit_type1, hit_type2;

  hit_type1 = m_CadInterface->shootRay(x, n, x_hit1, n_hit1, r_hit1);
  if (hit_type1 == CadInterface::Miss && !strict_direction) {
    n *= -1;
    hit_type1 = m_CadInterface->shootRay(x, n, x_hit1, n_hit1, r_hit1);
  }
  if (hit_type1 == CadInterface::Miss) {
    m_Failed = true;
    return x;
  }
  m_Failed = false;
  n *= -1;
  x_proj = x_hit1;
  m_LastNormal = n_hit1;
  m_LastRadius = r_hit1;
  if (!strict_direction) {
    hit_type2 = m_CadInterface->shootRay(x_hit1, n, x_hit2, n_hit2, r_hit2);
    if (hit_type2 != CadInterface::Miss) {
      if ((x - x_hit2).abs() < (x - x_hit1).abs()) {
        x_proj = x_hit2;
        m_LastNormal = n_hit2;
        m_LastRadius = r_hit2;
      }
    }
  }

  if (id_node != -1 && allow_search) {
    EG_VTKDCN(vtkDoubleArray, cl, m_FGrid, "node_meshdensity_desired");
    double L = 0.5*cl->GetValue(id_node);
    if ((x - x_proj).abs() > L) {
      vec3_t x_old;
      m_FGrid->GetPoint(id_node, x_old.data());
      double w = 0.1;
      vec3_t x_corr = w*x_proj + (1-w)*x_old;
      m_FGrid->GetPoints()->SetPoint(id_node, x_corr.data());
      x_corr = project(x_corr, id_node, true, m_FPart.globalNormal(id_node), false, false);
      m_FGrid->GetPoints()->SetPoint(id_node, x_old.data());
      if ((x_corr - x_proj).abs() > L) {
        m_Failed = true;
        return x;
      }
    }
  }
  /*
  {
    vec3_t x_old;
    vec3_t xp = x_proj;
    double scal, w=0;
    int count = 0;
    do {
      if (count >= 10) {
        m_Failed = true;
        return x;
      }
      vec3_t n1 = m_FPart.globalNormal(id_node);
      m_FGrid->GetPoint(id_node, x_old.data());
      m_FGrid->GetPoints()->SetPoint(id_node, x_proj.data());
      vec3_t n2 = m_FPart.globalNormal(id_node);
      m_FGrid->GetPoints()->SetPoint(id_node, x_old.data());
      scal = n1*n2;
      if (scal < 0.5) {
        w = min(1.0, w + 0.1);
        x_proj = w*x + (1-w)*xp;
      }
      ++count;
    } while (scal < 0.5);
  }
  */

  return x_proj;
}

double SurfaceProjection::getRadius(vtkIdType id_node)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  m_ForceRay = true;
  project(x, id_node);
  m_ForceRay = false;
  return m_LastRadius;
}

