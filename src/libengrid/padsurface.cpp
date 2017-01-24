// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2016 enGits GmbH                                      +
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

#include "padsurface.h"
#include "meshpartition.h"
#include "geometrytools.h"

PadSurface::PadSurface()
{
  m_Relative = true;
  m_Distance = 1.0;
  m_NewBC = 99;
}

void PadSurface::operate()
{
  MeshPartition part(m_Grid);
  part.setBCs(m_BCs);
  QList<QVector<vec3_t> > quads;
  foreach (vtkIdType id_cell, part.getCells()) {
    vec3_t n = GeometryTools::cellNormal(m_Grid, id_cell);
    n.normalise();
    EG_GET_CELL(id_cell, m_Grid);
    QVector<vec3_t> x(num_pts + 1);
    for (int i = 0; i < num_pts; ++i) {
      m_Grid->GetPoint(pts[i], x[i].data());
    }
    x[num_pts] = x[0];
    for (int i = 0; i < num_pts; ++i) {
      if (part.c2cGG(id_cell, i) < 0) {
        vec3_t b = x[i+1] - x[i];
        double L = m_Distance;
        if (m_Relative) {
          L = b.abs();
        }
        b.normalise();
        vec3_t a = b.cross(n);
        QVector<vec3_t> xn(4);
        xn[0] = x[i]   - L*b;
        xn[1] = x[i]   - L*b + L*a;
        xn[2] = x[i+1] + L*b + L*a;
        xn[3] = x[i+1] + L*b;
        quads << xn;
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, quads.size(), 4*quads.size());
  vtkIdType start_id = 0;
  vtkIdType pts[4];
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  foreach (QVector<vec3_t> quad, quads) {
    for (int i = 0; i < 4; ++i) {
      grid->GetPoints()->SetPoint(start_id + i, quad[i].data());
      pts[i] = start_id + i;
    }
    vtkIdType id_cell = grid->InsertNextCell(VTK_QUAD, 4, pts);
    cell_code->SetValue(id_cell, m_NewBC);
    start_id += 4;
  }
  MeshPartition new_part(grid, true);
  part.concatenatePartition(new_part);
}
