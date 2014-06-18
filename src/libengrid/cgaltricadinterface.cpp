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

#include "cgaltricadinterface.h"

CgalTriCadInterface::CgalTriCadInterface(vtkUnstructuredGrid *grid)
{
  m_Triangles.clear();
  QVector<vtkIdType> tris;
  getAllCellsOfType(VTK_TRIANGLE, tris, grid);
  m_Triangles.fill(Triangle(), tris.size());
  int i = 0;
  foreach (vtkIdType id_cell, tris) {
    EG_GET_CELL(id_cell, grid);
    vec3_t a, b, c;
    grid->GetPoint(pts[0], a.data());
    grid->GetPoint(pts[1], b.data());
    grid->GetPoint(pts[2], c.data());
    m_Triangles[i] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
    ++i;
  }
  m_Tree.rebuild(m_Triangles.begin(), m_Triangles.end());
  m_Tree.accelerate_distance_queries();
  setName("CgalTriCadInterface");
}

vec3_t CgalTriCadInterface::snap(vec3_t x, bool)
{
  Point p(x[0], x[1], x[2]);
  Point cp = m_Tree.closest_point(p);
  return vec3_t(cp[0], cp[1], cp[2]);
}
