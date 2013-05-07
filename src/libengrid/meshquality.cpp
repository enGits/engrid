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
#include "meshquality.h"

MeshQuality::MeshQuality()
{
}

void MeshQuality::computeNodesFromCells()
{
  EG_VTKDCN(vtkDoubleArray, node_mesh_quality, m_Grid, "node_mesh_quality");
  EG_VTKDCC(vtkDoubleArray, cell_mesh_quality, m_Grid, "cell_mesh_quality");
  EG_FORALL_NODES(id_node, m_Grid) {
    double mq = 1;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      mq = min(cell_mesh_quality->GetValue(id_cell), mq);
    }
    node_mesh_quality->SetValue(id_node, mq);
  }
}
