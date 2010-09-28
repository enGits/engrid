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
#include "deletetetras.h"

void DeleteTetras::operate()
{
  cout << "deleting tetrahedral cells" << endl;
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  QVector<vtkIdType> tetras, cells, nodes;

  {
    int N = 0;
    foreach (vtkIdType id_cell, m_Part.getCells()) {
      if (m_Grid->GetCellType(id_cell) == VTK_TETRA) {
        ++N;
      }
    }
    tetras.resize(N);
    N = 0;
    foreach (vtkIdType id_cell, m_Part.getCells()) {
      if (m_Grid->GetCellType(id_cell) == VTK_TETRA) {
        tetras[N] = id_cell;
        ++N;
      }
    }
  }

  getRestCells(m_Grid, tetras, cells);
  getNodesFromCells(cells, nodes, m_Grid);
  allocateGrid(new_grid, cells.size(), nodes.size());
  QVector<vtkIdType> old2new(m_Grid->GetNumberOfPoints(), -1);
  {
    vtkIdType id_new = 0;
    foreach (vtkIdType id_node, nodes) {  
      vec3_t x;
      m_Grid->GetPoints()->GetPoint(id_node, x.data());
      new_grid ->GetPoints()->SetPoint(id_new, x.data());
      copyNodeData(m_Grid, id_node, new_grid, id_new);
      old2new[id_node] = id_new;
      ++id_new;
    }
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType *pts, N_pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vtkIdType> new_pts(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      new_pts[i] = old2new[pts[i]];
    }
    vtkIdType cellType = m_Grid->GetCellType(id_cell);
    vtkIdType id_new = new_grid->InsertNextCell(cellType, N_pts, new_pts.data());
    copyCellData(m_Grid, id_cell, new_grid, id_new);
  }
  makeCopy(new_grid, m_Grid);
}

