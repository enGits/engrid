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

#include <QButtonGroup>
#include <QRadioButton>

///@@@  TODO: Change the checkboxes to a dropdown list or radiobuttons since the options are mutually exclusive anyway
void GuisetBoundaryCode::before()
{
  //read settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  ui.doubleSpinBoxFeatureAngle->setValue(local_qset.value("FeatureAngle", 45).toDouble());
  ui.spinBoxBoundaryCode->setValue(local_qset.value("BoundaryCode", 1).toInt());
  ui.checkBox_SelectAllVisible->setCheckState(int2CheckState(local_qset.value("SelectAllVisible", 0).toInt()));
  ui.checkBox_ProcessAll->setCheckState(int2CheckState(local_qset.value("ProcessAll", 0).toInt()));
  ui.checkBox_OnlyPickedCell->setCheckState(int2CheckState(local_qset.value("OnlyPickedCell", 0).toInt()));
  ui.checkBox_OnlyPickedCellAndNeighbours->setCheckState(int2CheckState(local_qset.value("OnlyPickedCellAndNeighbours", 0).toInt()));
  
  buttongroup = new QButtonGroup(this);
  radioButton_SelectAllVisible = new QRadioButton("Select all visible cells",this);
  radioButton_ProcessAll = new QRadioButton("Process all cells (even invisible ones)",this);
  radioButton_ProcessOnlyVisible = new QRadioButton("Process only visible cells",this);
  radioButton_OnlyPickedCell = new QRadioButton("Only picked cell",this);
  radioButton_OnlyPickedCellAndNeighbours = new QRadioButton("Only picked cell and neighbours",this);
  buttongroup->addButton(radioButton_SelectAllVisible,0);
  buttongroup->addButton(radioButton_ProcessAll,1);
  buttongroup->addButton(radioButton_ProcessOnlyVisible,2);
  buttongroup->addButton(radioButton_OnlyPickedCell,3);
  buttongroup->addButton(radioButton_OnlyPickedCellAndNeighbours,4);
  ui.verticalLayout_PickMethod->addWidget(radioButton_SelectAllVisible);
  ui.verticalLayout_PickMethod->addWidget(radioButton_ProcessAll);
  ui.verticalLayout_PickMethod->addWidget(radioButton_ProcessOnlyVisible);
  ui.verticalLayout_PickMethod->addWidget(radioButton_OnlyPickedCell);
  ui.verticalLayout_PickMethod->addWidget(radioButton_OnlyPickedCellAndNeighbours);
}

void GuisetBoundaryCode::operate()
{
  buttongroup->checkedId();
  cout<<"buttongroup->checkedId()="<<buttongroup->checkedId()<<endl;

  //save settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  local_qset.setValue("FeatureAngle", ui.doubleSpinBoxFeatureAngle->value());
  local_qset.setValue("BoundaryCode", ui.spinBoxBoundaryCode->value());
  local_qset.setValue("SelectAllVisible", ui.checkBox_SelectAllVisible->checkState());
  local_qset.setValue("ProcessAll", ui.checkBox_ProcessAll->checkState());
  local_qset.setValue("OnlyPickedCell", ui.checkBox_OnlyPickedCell->checkState());
  local_qset.setValue("OnlyPickedCellAndNeighbours", ui.checkBox_OnlyPickedCellAndNeighbours->checkState());
//   local_qset.setValue("PickType", ui.checkBox_OnlyPickedCellAndNeighbours->checkState());
  
  setBoundaryCode set_bc;
  if (GuiMainWindow::getPickedCell() >= 0) {
    
//     ui.radioButton_SelectAllVisible.
    
    set_bc.setFeatureAngle(ui.doubleSpinBoxFeatureAngle->value());
    set_bc.setBC(ui.spinBoxBoundaryCode->value());
    set_bc.setProcessAll(ui.checkBox_ProcessAll->checkState());
    set_bc.setSelectAllVisible(ui.checkBox_SelectAllVisible->checkState());
    set_bc.setOnlyPickedCell(ui.checkBox_OnlyPickedCell->checkState());
    set_bc.setOnlyPickedCellAndNeighbours(ui.checkBox_OnlyPickedCellAndNeighbours->checkState());
    
    cout<<"GuiMainWindow::getPickedCell()="<<GuiMainWindow::getPickedCell()<<endl;
    set_bc.setStart(GuiMainWindow::getPickedCell());
    set_bc();
  };
};
