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
#include "deletestraynodes.h"

DeleteStrayNodes::DeleteStrayNodes()
{
  EG_TYPENAME;
}

void DeleteStrayNodes::operate()
{
  cout << "removing 'stray' nodes" << endl;
  QVector<bool> active(m_Grid->GetNumberOfPoints(), false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      active[pts[i]] = true;
    }
  }
  int N = 0;
  foreach (bool is_active, active) {
    if (!is_active) {
      ++N;
    }
  }
  if (N == 1) {
    cout << "deleting 1 node!" << endl;
  } else {
    cout << "deleting " << N << " nodes!" << endl;
  }
  QVector<int> offset(m_Grid->GetNumberOfPoints(), 0);
  for (vtkIdType id_node = 1; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (!active[id_node]) {
      offset[id_node] = offset[id_node-1] + 1;
    } else {
      offset[id_node] = offset[id_node-1];
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints()-N);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (active[id_node]) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      new_grid->GetPoints()->SetPoint(id_node - offset[id_node], x.data());
      copyNodeData(m_Grid, id_node, new_grid, id_node - offset[id_node]);
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts, type;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    type = m_Grid->GetCellType(id_cell);
    QVector<vtkIdType> new_pts(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      new_pts[i] = pts[i] - offset[pts[i]];
    }
    vtkIdType id_new_cell = new_grid->InsertNextCell(type, N_pts, new_pts.data());
    copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
  }
  makeCopy(new_grid, m_Grid);
}

