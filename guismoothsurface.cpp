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
#include "guismoothsurface.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkLongArray.h>

void GuiSmoothSurface::before()
{
  populateBoundaryCodes(ui.listWidget, grid);
};

void GuiSmoothSurface::operate()
{
  QSet<int> bcs;
  getSelectedItems(ui.listWidget, bcs);
  QVector<vtkIdType> cells;
  getSurfaceCells(bcs, cells, grid);
  EG_VTKSP(vtkPolyData, pdata);
  addToPolyData(cells, pdata, grid);
  EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
  smooth->SetInput(pdata);
  smooth->FeatureEdgeSmoothingOn();
  smooth->SetFeatureAngle(ui.doubleSpinBox->value());
  smooth->SetNumberOfIterations(ui.spinBox->value());
  smooth->Update();
  EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
  for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
    vec3_t x;
    smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
    vtkIdType nodeId = node_index->GetValue(i);
    grid->GetPoints()->SetPoint(nodeId, x.data());
  };
};
