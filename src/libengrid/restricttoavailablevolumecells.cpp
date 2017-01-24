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
#include "restricttoavailablevolumecells.h"
#include "guimainwindow.h"

void RestrictToAvailableVolumeCells::operate()
{
  QList<vtkIdType> vol_cells, surf_cells;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      vol_cells.append(id_cell);
    }
  }
  foreach (vtkIdType id_cell, vol_cells) {
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
      if (id_neigh >= 0) {
        if (isSurface(id_neigh, m_Grid)) {
          surf_cells.append(id_neigh);
        }
      }
    }
  }
  MeshPartition vol_part(m_Grid);
  vol_part.setCells(vol_cells + surf_cells);
  EG_VTKSP(vtkUnstructuredGrid, new_grid1);
  vol_part.extractToVtkGrid(new_grid1);
  MeshPartition new_part1(new_grid1, true);
  QList<QVector<vtkIdType> > new_faces;
  for (vtkIdType id_cell = 0; id_cell < new_grid1->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, new_grid1)) {
      for (int i = 0; i < new_part1.c2cGSize(id_cell); ++i) {
        vtkIdType id_neigh = new_part1.c2cGG(id_cell, i);
        if (id_neigh == -1) {
          QVector<vtkIdType> pts;
          getFaceOfCell(new_grid1, id_cell, i, pts);
          new_faces.append(pts);
        }
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid2);
  allocateGrid(new_grid2, new_grid1->GetNumberOfCells() + new_faces.size(), new_grid1->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code1, new_grid1, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_code2, new_grid2, "cell_code");
  for (vtkIdType id_node = 0; id_node < new_grid1->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    new_grid1->GetPoint(id_node, x.data());
    new_grid2->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(new_grid1, id_node, new_grid2, id_node);
  }
  int bc_new = 0;
  for (vtkIdType id_cell = 0; id_cell < new_grid1->GetNumberOfCells(); ++id_cell) {
    copyCell(new_grid1, id_cell, new_grid2);
    copyCellData(new_grid1, id_cell, new_grid2, id_cell);
    if (isSurface(id_cell, new_grid1)) {
      bc_new = max(bc_new, cell_code1->GetValue(id_cell) + 1);
    } else {
      cell_code2->SetValue(id_cell, 0);
    }
  }
  foreach (QVector<vtkIdType> face, new_faces) {
    vtkIdType type = VTK_POLYGON;
    if (face.size() == 3) {
      type = VTK_TRIANGLE;
    } else if (face.size() == 4) {
      type = VTK_QUAD;
    }
    vtkIdType id_cell = new_grid2->InsertNextCell(type, face.size(), face.data());
    cell_code2->SetValue(id_cell, bc_new);
  }
  makeCopy(new_grid2, m_Grid);
  GuiMainWindow::pointer()->setBC(bc_new, BoundaryCondition("ami_blayer", "cyclic_AMI", bc_new));
  updateCellIndex(m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
}
