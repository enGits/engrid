// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#include "deletevolumegrid.h"

void DeleteVolumeGrid::operate()
{
  EG_VTKSP(vtkUnstructuredGrid, sgrid);
  QVector<vtkIdType> scells, snodes;
  QVector<int>       _snodes;
  getAllSurfaceCells(scells, m_Grid);
  getNodesFromCells(scells, snodes, m_Grid);
  createNodeMapping(snodes, _snodes, m_Grid);
  allocateGrid(sgrid, scells.size(), snodes.size());
  {
    vtkIdType newId = 0;
    foreach (vtkIdType nodeId, snodes) {  
      vec3_t x;
      m_Grid->GetPoint(nodeId, x.data());
      sgrid ->GetPoints()->SetPoint(newId, x.data());
      copyNodeData(m_Grid, nodeId, sgrid, newId);
      ++newId;
    }
  }
  foreach (vtkIdType cellId, scells) {
    vtkIdType *pts, Npts;
    m_Grid->GetCellPoints(cellId, Npts, pts);
    for (int i = 0; i < Npts; ++i) {
      pts[i] = _snodes[pts[i]];
    }
    vtkIdType cellType = m_Grid->GetCellType(cellId);
    vtkIdType newId = sgrid->InsertNextCell(cellType, Npts, pts);
    copyCellData(m_Grid, cellId, sgrid, newId);
  }
  makeCopy(sgrid, m_Grid);
}

