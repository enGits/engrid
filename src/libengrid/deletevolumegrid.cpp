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
#include "deletevolumegrid.h"

void DeleteVolumeGrid::operate()
{
  EG_VTKSP(vtkUnstructuredGrid, sgrid);
  QVector<vtkIdType> scells, snodes;
  getAllSurfaceCells(scells, m_Grid);
  getNodesFromCells(scells, snodes, m_Grid);
  allocateGrid(sgrid, scells.size(), snodes.size());
  {
    vtkIdType id_new = 0;
    foreach (vtkIdType id_node, snodes) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      sgrid ->GetPoints()->SetPoint(id_new, x.data());
      copyNodeData(m_Grid, id_node, sgrid, id_new);
      ++id_new;
    }
  }
  foreach (vtkIdType id_cell, scells) {
    vtkIdType id_new = copyCell(m_Grid, id_cell, sgrid);
    copyCellData(m_Grid, id_cell, sgrid, id_new);
  }
  makeCopy(sgrid, m_Grid);
}

