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
  //prepare radiobuttons
  m_ButtonGroup = new QButtonGroup(this);
  m_RadioButtonProcessOnlyVisible = new QRadioButton("Process only visible cells",this);
  m_RadioButtonProcessAll = new QRadioButton("Process all cells (even invisible ones)",this);
  m_RadioButtonSelectAllVisible = new QRadioButton("Select all visible cells",this);
  m_RadioButtonOnlyPickedCell = new QRadioButton("Only picked cell",this);
  m_RadioButtonOnlyPickedCellAndNeighbours = new QRadioButton("Only picked cell and neighbours",this);
  m_ButtonGroup->addButton(m_RadioButtonProcessOnlyVisible,0);
  m_ButtonGroup->addButton(m_RadioButtonProcessAll,1);
  m_ButtonGroup->addButton(m_RadioButtonSelectAllVisible,2);
  m_ButtonGroup->addButton(m_RadioButtonOnlyPickedCell,3);
  m_ButtonGroup->addButton(m_RadioButtonOnlyPickedCellAndNeighbours,4);
  ui.verticalLayout_PickMethod->addWidget(m_RadioButtonProcessOnlyVisible);
  ui.verticalLayout_PickMethod->addWidget(m_RadioButtonProcessAll);
  ui.verticalLayout_PickMethod->addWidget(m_RadioButtonSelectAllVisible);
  ui.verticalLayout_PickMethod->addWidget(m_RadioButtonOnlyPickedCell);
  ui.verticalLayout_PickMethod->addWidget(m_RadioButtonOnlyPickedCellAndNeighbours);
  
  //read settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  ui.doubleSpinBoxFeatureAngle->setValue(local_qset.value("FeatureAngle", 45).toDouble());
  ui.spinBoxBoundaryCode->setValue(local_qset.value("BoundaryCode", 1).toInt());
  m_ButtonGroup->button(local_qset.value("PickMethod",0).toInt())->setChecked(true);
}

void GuiSetBoundaryCode::operate()
{
  m_ButtonGroup->checkedId();
  cout<<"buttongroup->checkedId()="<<m_ButtonGroup->checkedId()<<endl;

  //save settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  local_qset.setValue("FeatureAngle", ui.doubleSpinBoxFeatureAngle->value());
  local_qset.setValue("BoundaryCode", ui.spinBoxBoundaryCode->value());
  local_qset.setValue("PickMethod", m_ButtonGroup->checkedId());
  
  SetBoundaryCode set_bc;
  set_bc.setGrid(grid);
  set_bc.setAllSurfaceCells();
  if (mainWindow()->getPickedCell() >= 0) {
    set_bc.setFeatureAngle(ui.doubleSpinBoxFeatureAngle->value());
    set_bc.setBC(ui.spinBoxBoundaryCode->value());
    
    set_bc.setProcessAll(m_ButtonGroup->button(1)->isChecked());
    set_bc.setSelectAllVisible(m_ButtonGroup->button(2)->isChecked());
    set_bc.setOnlyPickedCell(m_ButtonGroup->button(3)->isChecked());
    set_bc.setOnlyPickedCellAndNeighbours(m_ButtonGroup->button(4)->isChecked());
    
    cout << "GuiMainWindow::getPickedCell()=" << mainWindow()->getPickedCell() << endl;
    set_bc.setStart(mainWindow()->getPickedCell());
    set_bc();
  }
}
