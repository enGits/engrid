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
#include "meshquality.h"

MeshQuality::MeshQuality()
{
  m_Name = "unknown mesh quality";
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

void MeshQuality::printCellInfo(int indent)
{
  double min_quality = 1e99;
  double max_quality = 0;
  double mean_quality = 0;
  EG_VTKDCC(vtkDoubleArray, cell_mesh_quality, m_Grid, "cell_mesh_quality");
  EG_FORALL_CELLS(id_cell, m_Grid) {
    double q = cell_mesh_quality->GetValue(id_cell);
    min_quality = min(q, min_quality);
    max_quality = max(q, max_quality);
    mean_quality += q;
  }
  mean_quality /= m_Grid->GetNumberOfCells();
  QString indent_str;
  indent_str.fill(' ', indent);
  cout << qPrintable(indent_str) << qPrintable(name()) << endl;
  cout << qPrintable(indent_str) << "min  : " << min_quality << endl;
  cout << qPrintable(indent_str) << "mean : " << mean_quality << endl;
  cout << qPrintable(indent_str) << "max  : " << max_quality << endl;
}
