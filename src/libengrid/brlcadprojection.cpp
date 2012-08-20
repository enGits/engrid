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
  EG_VTKDCN(vtkDoubleArray, cl, m_FGrid, "node_meshdensity_desired");
  m_Failed = true;
  if (id_node != -1) {
    setEpsilon(0.1*cl->GetValue(id_node));
  }
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
    cout << "vector defect (id_node=" << id_node << endl;
    return x;
    EG_BUG;
  }

  vec3_t x_proj = x;
  m_LastNormal = n;
  m_LastRadius = 1e10;

  vec3_t x_hit, n_hit;
  double r_hit;

  PositionType pos_type = position(x - 0.001*n, n);
  if (pos_type == Surface && !m_ForceRay) {
    //return x;
    if (id_node != -1) {
      x -= 0.01*cl->GetValue(id_node)*n;
    }
  }
  if (pos_type == Outside) {
    n *= -1;
  }

  // first shot along the provided (or computed) mesh normal
  HitType hit_type = shootRay(x - 0.01*n, n, x_hit, n_hit, r_hit);
  if (hit_type != BrlCadInterface::Miss) {
    x_proj = x_hit;
    m_LastNormal = n_hit;
    m_LastRadius = r_hit;
    m_Failed = false;
  }

  /*
  if (x_proj.abs() > 0.99) {
    if (fabs(x_proj.abs() - 1) > 1e-4) {
      cout << x_proj.abs() << " " << x.abs() << " " << pos_type << "," << hit_type << endl;
    }
  }
  */

  setEpsilon(0);
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
