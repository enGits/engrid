//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "guimirrormesh.h"

void GuiMirrorMesh::operate()
{
  vec3_t x0, n;
  x0[0] = ui.lineEditX0->text().toDouble();
  x0[1] = ui.lineEditY0->text().toDouble();
  x0[2] = ui.lineEditZ0->text().toDouble();
  n[0] = ui.lineEditNX->text().toDouble();
  n[1] = ui.lineEditNY->text().toDouble();
  n[2] = ui.lineEditNZ->text().toDouble();
  n.normalise();
  EG_VTKSP(vtkUnstructuredGrid, mirror_grid);
  EgVtkObject::allocateGrid(mirror_grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    double h = n*(x-x0);
    vec3_t xm = x - 2*h*n;
    mirror_grid->GetPoints()->SetPoint(id_node, xm.data());
    copyNodeData(m_Grid, id_node, mirror_grid, id_node);
  }
  vtkIdType id_new_cell;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType new_pts[N_pts];
    vtkIdType cell_type = m_Grid->GetCellType(id_cell);
    if (cell_type == VTK_WEDGE) {
      new_pts[0] = pts[3];
      new_pts[1] = pts[4];
      new_pts[2] = pts[5];
      new_pts[3] = pts[0];
      new_pts[4] = pts[1];
      new_pts[5] = pts[2];
    } else if (cell_type == VTK_HEXAHEDRON) {
      new_pts[0] = pts[4];
      new_pts[1] = pts[5];
      new_pts[2] = pts[6];
      new_pts[3] = pts[7];
      new_pts[4] = pts[0];
      new_pts[5] = pts[1];
      new_pts[6] = pts[2];
      new_pts[7] = pts[3];
    } else if (cell_type == VTK_PYRAMID) {
      EG_BUG;
    } else {
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = pts[N_pts - i - 1];
      }
    }
    id_new_cell = mirror_grid->InsertNextCell(m_Grid->GetCellType(id_cell), N_pts, new_pts);
    copyCellData(m_Grid, id_cell, mirror_grid, id_new_cell);
  }
  MeshPartition part1(m_Grid, true);
  MeshPartition part2(mirror_grid, true);
  part1.addPartition(part2);
  eliminateDuplicateCells();
}

