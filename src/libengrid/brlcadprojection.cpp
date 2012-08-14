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
    cout << "vector defect (id_node=" << id_node << endl;
    return x;
    EG_BUG;
  }

  vec3_t x_proj = x;

  vec3_t x_hit, n_hit;
  double r_hit;

  PositionType pos_type = position(x, n);
  if (pos_type == Surface && !m_ForceRay) {
    //return x;
  }
  if (pos_type == Outside) {
    n *= -1;
  }

  // first shot along the provided (or computed) mesh normal
  if (shootRay(x, n, x_hit, n_hit, r_hit) != BrlCadInterface::Miss) {
    x_proj = x_hit;
    m_LastNormal = n_hit;
    m_LastRadius = r_hit;

    /*
    // iterate into corners
    int iteration_count = 0;
    double L_max = (x - x_proj).abs();
    n = n_hit;
    vec3_t x0 = x;
    do {
      if (shootRay(x0, n, x_hit, n_hit, r_hit) != BrlCadInterface::Miss) {
        double L = (x_hit - x).abs();
        if (L < L_max) {
          x_proj = x_hit;
          m_LastNormal = n_hit;
          m_LastRadius = r_hit;
          break;
        }
      }
      x0 = 0.5*(x0 + x_proj);
      ++iteration_count;
    } while (iteration_count < 10);
    */
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
