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
#include <vtkMergePoints.h>

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

void MeshPartition::extractToVtkGrid(vtkUnstructuredGrid *new_grid)
{
  allocateGrid(new_grid, cells.size(), nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vec3_t x;
    grid->GetPoints()->GetPoint(nodes[i_nodes], x.data());
    new_grid->GetPoints()->SetPoint(i_nodes, x.data());
    copyNodeData(grid, nodes[i_nodes], new_grid, i_nodes);
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = grid->GetCellType(id_cell);
    grid->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, pts);
    copyCellData(grid, id_cell, new_grid, id_new_cell);
  }
}

void MeshPartition::addPartition(const MeshPartition& part)
{
  if (grid == part.grid) {
  } else {
    EG_VTKSP(vtkUnstructuredGrid, new_grid);
    EG_VTKSP(vtkMergePoints,loc);
    loc->SetDataSet(grid);
    loc->BuildLocator();

    QVector<vtkIdType> pnode2node(part.grid->GetNumberOfPoints());
    vtkIdType N = grid->GetNumberOfPoints();
    foreach (vtkIdType id_pnode, part.nodes) {
      vec3_t x;
      part.grid->GetPoint(id_pnode, x.data());
      if (loc->IsInsertedPoint(x.data())) {
        pnode2node[id_pnode] = loc->FindClosestPoint(x.data());
      } else {
        pnode2node[id_pnode] = N;
        ++N;
      }
    }
    allocateGrid(new_grid, grid->GetNumberOfCells() + part.cells.size(), N);
    for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      grid->GetPoint(id_node, x.data());
      new_grid->GetPoints()->SetPoint(id_node, x.data());
      copyNodeData(grid, id_node, new_grid, id_node);
    }
    foreach (vtkIdType id_pnode, part.nodes) {
      vec3_t x;
      part.grid->GetPoint(id_pnode, x.data());
      new_grid->GetPoints()->SetPoint(pnode2node[id_pnode], x.data());
      copyNodeData(part.grid, id_pnode, new_grid, pnode2node[id_pnode]);
    }
    QList<vtkIdType> new_cells;
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = grid->GetCellType(id_cell);
      grid->GetCellPoints(id_cell, N_pts, pts);
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, pts);
      copyCellData(grid, id_cell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    foreach (vtkIdType id_pcell, part.cells) {
      vtkIdType N_pts, *pts;
      vtkIdType new_pts[N_pts];
      vtkIdType type_cell = part.grid->GetCellType(id_pcell);
      part.grid->GetCellPoints(id_pcell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = pnode2node[pts[i]];
      }
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, new_pts);
      copyCellData(part.grid, id_pcell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    makeCopy(new_grid, grid);
  }
}
