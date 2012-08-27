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

#include "brlcadprojection.h"

BrlCadProjection::BrlCadProjection(QString file_name, QString object_name)
{
  setupBrlCad(file_name, object_name);
  m_ForceRay = false;
  m_Failed = false;
}

BrlCadProjection::~BrlCadProjection()
{

}

vec3_t BrlCadProjection::project(vec3_t x, vtkIdType id_node, bool, vec3_t v)
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
  HitType hit_type1, hit_type2;

  hit_type1 = shootRay(x, n, x_hit1, n_hit1, r_hit1);
  if (hit_type1 == Miss) {
    n *= -1;
    hit_type1 = shootRay(x, n, x_hit1, n_hit1, r_hit1);
  }
  if (hit_type1 == Miss) {
    m_Failed = true;
    return x;
  }
  m_Failed = false;
  n *= -1;
  hit_type2 = shootRay(x_hit1, n, x_hit2, n_hit2, r_hit2);
  x_proj = x_hit1;
  m_LastNormal = n_hit1;
  m_LastRadius = r_hit1;
  if (hit_type2 != Miss) {
    if ((x - x_hit2).abs() < (x - x_hit1).abs()) {
      x_proj = x_hit2;
      m_LastNormal = n_hit2;
      m_LastRadius = r_hit2;
    }
  }
  return x_proj;
}

double BrlCadProjection::getRadius(vtkIdType id_node)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  m_ForceRay = true;
  project(x, id_node);
  m_ForceRay = false;
  return m_LastRadius;
}
