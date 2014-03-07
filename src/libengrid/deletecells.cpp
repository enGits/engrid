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
#include "deletecells.h"

void DeleteCells::setTraceCells(const QVector<vtkIdType> &cells)
{
  trace_cells.resize(cells.size()); 
  qCopy(cells.begin(), cells.end(), trace_cells.begin()); 
}

void DeleteCells::getTraceCells(QVector<vtkIdType> &cells)
{
  cells.resize(trace_cells.size()); 
  qCopy(trace_cells.begin(), trace_cells.end(), cells.begin()); 
}

void DeleteCells::operate()
{
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  QVector<vtkIdType> new_cells;
  QVector<vtkIdType> new_nodes;
  getRestCells(m_Grid, del_cells, new_cells);
  getNodesFromCells(new_cells, new_nodes, m_Grid);
  allocateGrid(new_grid, new_cells.size(), new_nodes.size());
  QVector<vtkIdType> old2new_nodes(m_Grid->GetNumberOfPoints(), -1);
  QVector<vtkIdType> old2new_cells(m_Grid->GetNumberOfCells(), -1);
  {
    vtkIdType id_new = 0;
    foreach (vtkIdType id_node, new_nodes) {
      vec3_t x;
      m_Grid->GetPoints()->GetPoint(id_node, x.data());
      new_grid ->GetPoints()->SetPoint(id_new, x.data());
      copyNodeData(m_Grid, id_node, new_grid, id_new);
      old2new_nodes[id_node] = id_new;
      ++id_new;
    }
  }
  {
    foreach (vtkIdType id_cell, new_cells) {
      vtkIdType *pts, N_pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      QVector<vtkIdType> new_pts(N_pts);
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = old2new_nodes[pts[i]];
      }
      vtkIdType cellType = m_Grid->GetCellType(id_cell);
      vtkIdType id_new = new_grid->InsertNextCell(cellType, N_pts, new_pts.data());
      copyCellData(m_Grid, id_cell, new_grid, id_new);
      old2new_cells[id_cell] = id_new;
    }
    QList<vtkIdType> new_trace_cells;
    foreach (vtkIdType id_cell, trace_cells) {
      if (old2new_cells[id_cell] != -1) {
        new_trace_cells.append(old2new_cells[id_cell]);
      }
    }
    trace_cells.resize(new_trace_cells.size());
    qCopy(new_trace_cells.begin(), new_trace_cells.end(), trace_cells.begin());
  }
  makeCopy(new_grid, m_Grid);
}

