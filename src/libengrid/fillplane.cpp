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

#include "fillplane.h"

vec3_t FillPlane::toPlane(vec3_t x)
{
  x -= m_X0;
  x -= (x*m_N)*m_N;
  return vec3_t(x*m_G1, x*m_G2, 0);
}

vec3_t FillPlane::fromPlane(vec3_t x)
{
  return m_X0 + x[0]*m_G1 + x[1]*m_G2;
}

void FillPlane::createEdgesOnPlane(vtkUnstructuredGrid *edge_grid)
{

}

void FillPlane::operate()
{
  m_G1 = GeometryTools::orthogonalVector(m_N);
  m_G2 = m_N.cross(m_G1);
  m_N.normalise();
  m_G1.normalise();
  m_G2.normalise();
  EG_VTKSP(vtkUnstructuredGrid, edge_grid);
  createEdgesOnPlane(edge_grid);
}


