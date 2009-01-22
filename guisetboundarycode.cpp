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
#include "guisetboundarycode.h"
#include "setboundarycode.h"
#include "guimainwindow.h"

void GuiSetBoundaryCode::operate()
{
  SetBoundaryCode set_bc;
  if (GuiMainWindow::getPickedCell() >= 0) {
    set_bc.setFeatureAngle(ui.doubleSpinBoxFeatureAngle->value());
    set_bc.setBC(ui.spinBoxBoundaryCode->value());
    set_bc.setStart(GuiMainWindow::getPickedCell());
    set_bc();
  };
};
