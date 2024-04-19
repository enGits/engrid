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
#include "correctsurfaceorientation.h"
#include <vtkSmartPointer.h>

void CorrectSurfaceOrientation::pass1()
{
  for (int i = 0; i < pair.size(); ++i) {
    pair[i].terminate = false;
    if (pair[i].item1 >= 0) {
      vtkSmartPointer<vtkIdList> pts1 = vtkSmartPointer<vtkIdList>::New();
      vtkSmartPointer<vtkIdList> pts2 = vtkSmartPointer<vtkIdList>::New();
      m_Grid->GetCellPoints(pair[i].item1, pts1);
      m_Grid->GetCellPoints(pair[i].item2, pts2);
      vtkIdType num_pts1 = pts1->GetNumberOfIds();
      vtkIdType num_pts2 = pts2->GetNumberOfIds();

      bool ok = false;
      for (int j1 = 0; j1 < num_pts1; ++j1) {
        for (int j2 = 0; j2 < num_pts2; ++j2) {
          vtkIdType node_11 = pts1->GetId(j1);
          vtkIdType node_21 = pts2->GetId(j2);
          vtkIdType node_12, node_22;
          if (j1 < num_pts1 - 1) node_12 = pts1->GetId(j1 + 1);
          else                   node_12 = pts1->GetId(0);
          if (j2 < num_pts2 - 1) node_22 = pts2->GetId(j2 + 1);
          else                   node_22 = pts2->GetId(0);
          if ((node_11 == node_22) && (node_12 == node_21)) {
            ok = true;
            break;
          }
        }
      }
      if (!ok) {
        m_Grid->GetCells()->ReverseCellAtId(pair[i].item2);
        // QVector<vtkIdType> nodes(num_pts2);
        // for (vtkIdType j = 0; j < num_pts2; ++j) nodes[j] = pts2->GetId(j);
        // for (vtkIdType j = 0; j < num_pts2; ++j) pts2->SetId(num_pts2 - j - 1, nodes[j]);
      }
    }
  }
}

