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

void GuiSetBoundaryCode::before()
{
  //read settings
  QSettings local_qset("enGits","enGrid_GuiSetBoundaryCode");
  ui.doubleSpinBoxFeatureAngle->setValue(local_qset.value("FeatureAngle", 45).toDouble());
  ui.spinBoxBoundaryCode->setValue(local_qset.value("BoundaryCode", 1).toInt());
  ui.checkBox_SelectAllVisible->setCheckState(int2CheckState(local_qset.value("SelectAllVisible", 0).toInt()));
  ui.checkBox_ProcessAll->setCheckState(int2CheckState(local_qset.value("ProcessAll", 0).toInt()));
}

void GuiSetBoundaryCode::operate()
{
  //save settings
  QSettings local_qset("enGits","enGrid_GuiSetBoundaryCode");
  local_qset.setValue("FeatureAngle", ui.doubleSpinBoxFeatureAngle->value());
  local_qset.setValue("BoundaryCode", ui.spinBoxBoundaryCode->value());
  local_qset.setValue("SelectAllVisible", ui.checkBox_SelectAllVisible->checkState());
  local_qset.setValue("ProcessAll", ui.checkBox_ProcessAll->checkState());
  
  SetBoundaryCode set_bc;
  if (GuiMainWindow::getPickedCell() >= 0) {
    set_bc.setFeatureAngle(ui.doubleSpinBoxFeatureAngle->value());
    set_bc.setBC(ui.spinBoxBoundaryCode->value());
    set_bc.setProcessAll(ui.checkBox_ProcessAll->checkState());
    set_bc.setSelectAllVisible(ui.checkBox_SelectAllVisible->checkState());
    set_bc.setStart(GuiMainWindow::getPickedCell());
    set_bc();
  };
};
