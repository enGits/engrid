//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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

#include "meshpartition.h"

MeshPartition::MeshPartition()
{
  grid = NULL;
  original_orientation = true;
}

MeshPartition::MeshPartition(QString volume_name)
{
  grid = GuiMainWindow::pointer()->getGrid();
  VolumeDefinition V = GuiMainWindow::pointer()->getVol(volume_name);
  QList<vtkIdType> cls;
  QList<vtkIdType> reo_fcs;
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell)
  {
    if (isSurface(id_cell, grid)) {
      int bc = cell_code->GetValue(id_cell);
      if (V.getSign(bc) != 0) {
        cls.append(id_cell);
        if (V.getSign(bc) == -1) {
          reo_fcs.append(id_cell);
        }
      }
    } else {
      if (cell_code->GetValue(id_cell) == V.getVC()) {
        cls.append(id_cell);
      }
    }
  }
  re_orientate_faces.resize(reo_fcs.size());
  qCopy(reo_fcs.begin(), reo_fcs.end(), re_orientate_faces.begin());
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  createCellMapping(cells, _cells, grid);
  createNodeMapping(nodes, _nodes, grid);
  original_orientation = true;
}

void MeshPartition::changeOrientation()
{
  foreach (vtkIdType id_cell, re_orientate_faces) {
    vtkIdType Npts, *pts;
    grid->GetCellPoints(id_cell, Npts, pts);
    vtkIdType new_pts[Npts];
    for (int i = 0; i < Npts; ++i) {
      new_pts[i] = pts[Npts-1-i];
    }
    grid->ReplaceCell(id_cell, Npts, new_pts);
  }
  original_orientation = !original_orientation;
}

void MeshPartition::setVolumeOrientation()
{
  if (original_orientation) {
    changeOrientation();
  }
}

void MeshPartition::setOriginalOrientation()
{
  if (!original_orientation) {
    changeOrientation();
  }
}

void MeshPartition::setRemainder(const MeshPartition& part)
{
  setGrid(part.getGrid());
  QVector<vtkIdType> rcells;
  getRestCells(grid, cells, rcells);
  setCells(rcells);
}

