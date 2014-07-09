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
#include "checksurfaceintegrity.h"
// #include "egvtkobject.h"
#include "geometrytools.h"
using namespace GeometryTools;

CheckSurfaceIntegrity::CheckSurfaceIntegrity() : SurfaceOperation()
{
  EG_TYPENAME;
}

void CheckSurfaceIntegrity::operate()
{
  bool ok = isWaterTight();
  cout << endl;
  if (ok) {
    cout << "The mesh is OK!" << endl;
  } else {
    cout << "The mesh has defects!" << endl;
    for (int i = 0; i < m_NumCells.size(); ++i) {
      if (m_NumCells[i] > 0) {
        cout << "  " << m_NumCells[i] << " edges with " << i << " faces" << endl;
      }
    }
    cout << "  " << m_NumEmptyCells << " faces are degenerated" << endl;
    cout << endl;
  }
}

bool CheckSurfaceIntegrity::isWaterTight()
{
  setAllSurfaceCells();
  m_BadCells.clear();
  m_NumCells.fill(0, 1000);

  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();
  l2g_t  cells = getPartCells();
  
  bool first = true;
  foreach(vtkIdType id_node1, nodes) {
    vec3_t x;
    m_Grid->GetPoints()->GetPoint(id_node1, x.data());
    if(!checkVector(x)) {
      qDebug()<<"point "<<id_node1<<" is NAN or INF";
    };
    foreach(int i_node2, n2n[_nodes[id_node1]]) {
      vtkIdType id_node2 = nodes[i_node2];
      QSet <vtkIdType> edge_cells;
      int N = getEdgeCells(id_node1, id_node2, edge_cells);
      if(first) {
        first = false;
        m_NumMin = N;
        m_NumMax = N;
      } else {
        m_NumMin = min(m_NumMin, N);
        m_NumMax = max(m_NumMax, N);
      }
      if (N >= 1000) {
        EG_BUG;
      }
      ++m_NumCells[N];
      if (edge_cells.size() != 2) {
        m_BadCells.unite(edge_cells);
      }
    }
  }
  
  m_NumEmptyCells = 0;
  foreach(vtkIdType id_cell, cells) {
    if (cellVA(m_Grid, id_cell) <= 0) {
      ++m_NumEmptyCells;
    }
  }
  
  if (m_NumMin == 2 && m_NumMax == 2 && m_NumEmptyCells == 0) {
    return(true);
  } else {
    return(false);
  }
}
