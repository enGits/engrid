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

#include "seligairfoilreader.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"

#include <vtkDelaunay2D.h>

SeligAirfoilReader::SeligAirfoilReader()
{
  setFormat("Selig airfoil data file (*.dat *.DAT)");
}

void SeligAirfoilReader::operate()
{
  try {
    readInputFileName("");
    if (isValid()) {
      QFile file(getFileName());
      file.open(QIODevice::ReadOnly | QIODevice::Text);
      QTextStream f(&file);
      f.readLine();
      EG_VTKSP(vtkPolyData, poly);
      double read_value;
      int num_upper, num_lower;
      f >> read_value;
      num_upper = int(read_value);
      f >> read_value;
      num_lower = int(read_value);
      int num_nodes = num_lower + num_upper - 1;
      QVector<vec3_t> coord(num_nodes, vec3_t(0,0,0));
      for (int i = 0; i < num_upper; ++i) {
        f >> coord[i][0] >> coord[i][1];
      }
      double dummy;
      f >> dummy;
      f >> dummy;
      for (int i = num_nodes - 1; i >= num_upper; --i) {
        f >> coord[i][0] >> coord[i][1];
      }
      EG_VTKSP(vtkPoints, points);
      points->SetNumberOfPoints(num_nodes);
      poly->SetPoints(points);
      poly->Allocate(num_nodes);
      for (vtkIdType id_node = 0; id_node < num_nodes; ++id_node) {
        poly->GetPoints()->SetPoint(id_node, coord[id_node].data());
      }
      for (vtkIdType id_node = 0; id_node < num_nodes; ++id_node) {
        vtkIdType pts[2];
        pts[0] = id_node;
        if (id_node < num_nodes-1) {
          pts[1] = id_node + 1;
        } else {
          pts[1] = 0;
        }
        poly->InsertNextCell(VTK_LINE, 2, pts);
      }

      EG_VTKSP(vtkDelaunay2D, tri);
      tri->SetInputData(poly);
      tri->SetSourceData(poly);
      EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, poly2ug);
      poly2ug->SetInputConnection(tri->GetOutputPort());
      poly2ug->Update();
      makeCopy(poly2ug->GetOutput(), m_Grid);
      updateNodeIndex(m_Grid);
      updateCellIndex(m_Grid);
    }
  } catch (Error err) {
    err.display();
  }
}

