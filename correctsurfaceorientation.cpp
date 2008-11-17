//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "correctsurfaceorientation.h"

void CorrectSurfaceOrientation::pass1()
{
  for (int i = 0; i < pair.size(); ++i) {
    pair[i].terminate = false;
    if (pair[i].item1 >= 0) {
      vtkIdType Npts1, *pts1;
      vtkIdType Npts2, *pts2;
      grid->GetCellPoints(pair[i].item1, Npts1, pts1);
      grid->GetCellPoints(pair[i].item2, Npts2, pts2);
      bool ok = false;
      for (int j1 = 0; j1 < Npts1; ++j1) {
        for (int j2 = 0; j2 < Npts2; ++j2) {
          vtkIdType node_11 = pts1[j1];
          vtkIdType node_21 = pts2[j2];
          vtkIdType node_12, node_22;
          if (j1 < Npts1 - 1) node_12 = pts1[j1 + 1];
          else                node_12 = pts1[0];
          if (j2 < Npts2 - 1) node_22 = pts2[j2 + 1];
          else                node_22 = pts2[0];
          if ((node_11 == node_22) && (node_12 == node_21)) {
            ok = true;
            break;
          };
        };
      };
      if (!ok) {
        QVector<vtkIdType> nodes(Npts2);
        for (vtkIdType j = 0; j < Npts2; ++j) nodes[j]            = pts2[j];
        for (vtkIdType j = 0; j < Npts2; ++j) pts2[Npts2 - j - 1] = nodes[j];
      };
    };
  };
};

