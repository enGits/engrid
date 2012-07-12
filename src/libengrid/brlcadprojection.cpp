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
  m_Rtip = rt_dirbuild(qPrintable(file_name), m_IdBuf, sizeof(m_IdBuf));
  if (m_Rtip == RTI_NULL) {
    EG_ERR_RETURN("Unable to open BRL-CAD database!");
  }
  if (rt_gettree(m_Rtip, qPrintable(object_name)) < 0) {
    EG_ERR_RETURN("unable to access object \"" + object_name + "\"");
  }

  application ap = {0};
  m_Ap = ap;
  m_Ap.a_rt_i   = m_Rtip;
  m_Ap.a_onehit = 0; // no X-ray functionality

  rt_prep_parallel(m_Rtip, 1);
}

BrlCadProjection::~BrlCadProjection()
{

}

vec3_t BrlCadProjection::project(vec3_t x, vtkIdType id_node)
{
  vec3_t n = m_FPart.globalNormal(id_node);
  /*
  if (!checkVector(x)) {
    EG_BUG;
  }
  bool first = true;
  VSET(m_Ap.a_ray.r_pt,  x[0], x[1], x[2]);
  VSET(m_Ap.a_ray.r_dir, v[0], v[1], v[2]);
  m_Hit = false;
  m_Ap.a_hit  = CreateBrlCadTesselation::hit;
  m_Ap.a_miss = CreateBrlCadTesselation::miss;
  rt_shootray(&m_Ap);
  x_in = m_XIn;
  n_in = m_InNormal;
  x_out = m_XOut;
  n_out = m_OutNormal;
  x = x_out;
  if (!checkVector(x_in)) {
    EG_BUG;
  }
  if (!checkVector(x_out)) {
    EG_BUG;
  }
  if (!checkVector(n_in)) {
    EG_BUG;
  }
  if (!checkVector(n_out)) {
    EG_BUG;
  }
  return m_Hit;
  */
}

vec3_t BrlCadProjection::projectFree(vec3_t x, vtkIdType id_node, bool)
{
  project(x, id_node);
}

vec3_t BrlCadProjection::projectRestricted(vec3_t x, vtkIdType id_node, bool)
{
  project(x, id_node);
}
