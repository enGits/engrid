//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "guiimproveaspectratio.h"
#include "vtkEgEliminateShortEdges.h"

void GuiImproveAspectRatio::before()
{
  ui.lineEditAspectRatio->setText("100.0");
};

void GuiImproveAspectRatio::operate()
{
  EG_VTKSP(vtkEgEliminateShortEdges, elem);
  int N_elim = 0;
  int N_sweeps = 0;
  elem->SetMaxRatio(ui.lineEditAspectRatio->text().toDouble());
  elem->SetMaxLength(ui.lineEditLength->text().toDouble());
  EG_VTKSP(vtkUnstructuredGrid,ug);
  do {
    ug->DeepCopy(grid);
    elem->SetInput(ug);
    elem->Update();
    grid->DeepCopy(elem->GetOutput());
    N_elim += elem->GetNumEliminated();
    ++N_sweeps;
  } while (elem->GetNumEliminated() > 0);
  QString text, t1, t2;
  t1.setNum(N_elim);
  t2.setNum(N_sweeps);
};

