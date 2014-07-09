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
#include "guisetboundarycode.h"
#include "setboundarycode.h"
#include "guimainwindow.h"

///\todo Add automatic BC generation using feature edges
void GuiSetBoundaryCode::before()
{
  //prepare radiobuttons
  m_ButtonGroup = new QButtonGroup(this);
  m_RadioButtonProcessOnlyVisible = new QRadioButton("Process only visible cells",this);
  m_RadioButtonProcessAll = new QRadioButton("Process all cells (even invisible ones)",this);
  m_RadioButtonSelectAllVisible = new QRadioButton("Select all visible cells",this);
  m_RadioButtonOnlyPickedCell = new QRadioButton("Only picked cell",this);
  m_RadioButtonOnlyPickedCellAndNeighbours = new QRadioButton("Only picked cell and neighbours",this);
  m_RadioButtonAuto = new QRadioButton("automatic",this);
  m_ButtonGroup->addButton(m_RadioButtonProcessOnlyVisible,0);
  m_ButtonGroup->addButton(m_RadioButtonProcessAll,1);
  m_ButtonGroup->addButton(m_RadioButtonSelectAllVisible,2);
  m_ButtonGroup->addButton(m_RadioButtonOnlyPickedCell,3);
  m_ButtonGroup->addButton(m_RadioButtonOnlyPickedCellAndNeighbours,4);
  m_ButtonGroup->addButton(m_RadioButtonAuto,5);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonProcessOnlyVisible);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonProcessAll);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonSelectAllVisible);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonOnlyPickedCell);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonOnlyPickedCellAndNeighbours);
  m_Ui.verticalLayout_PickMethod->addWidget(m_RadioButtonAuto);

  //read settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  m_Ui.doubleSpinBoxFeatureAngle->setValue(local_qset.value("FeatureAngle", 45).toDouble());
  m_Ui.spinBoxBoundaryCode->setValue(local_qset.value("BoundaryCode", 1).toInt());
  m_ButtonGroup->button(local_qset.value("PickMethod",0).toInt())->setChecked(true);
}

void GuiSetBoundaryCode::operate()
{
  m_ButtonGroup->checkedId();
  cout<<"buttongroup->checkedId()="<<m_ButtonGroup->checkedId()<<endl;

  //save settings
  QSettings local_qset("enGits","enGrid_GuisetBoundaryCode");
  local_qset.setValue("FeatureAngle", m_Ui.doubleSpinBoxFeatureAngle->value());
  local_qset.setValue("BoundaryCode", m_Ui.spinBoxBoundaryCode->value());
  local_qset.setValue("PickMethod", m_ButtonGroup->checkedId());
  
  SetBoundaryCode set_bc;
  set_bc.setGrid(m_Grid);
  set_bc.setAllSurfaceCells();
  if (m_RadioButtonAuto->isChecked()) {
    QSet <int> display_bcs;
    GuiMainWindow::pointer()->getDisplayBoundaryCodes(display_bcs);
    int bc = m_Ui.spinBoxBoundaryCode->value();
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      bc = max(bc, cell_code->GetValue(id_cell));
      if (display_bcs.contains(cell_code->GetValue(id_cell))) {
        cell_code->SetValue(id_cell, 9999);
      }
    }
    bool done = false;
    do {
      vtkIdType id_start = -1;
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        if (cell_code->GetValue(id_cell) == 9999) {
          id_start = id_cell;
          break;
        }
      }
      if (id_start == -1) {
        done = true;
      } else {
        set_bc.setFeatureAngle(m_Ui.doubleSpinBoxFeatureAngle->value());
        set_bc.setBC(bc);
        set_bc.setProcessAll(false);
        set_bc.setSelectAllVisible(false);
        set_bc.setOnlyPickedCell(false);
        set_bc.setOnlyPickedCellAndNeighbours(false);
        set_bc.setStart(id_start);
        QString bc_name = GuiMainWindow::pointer()->getBC(bc).getName();
        if (bc_name == "unknown") {
          bc_name.setNum(bc);
          bc_name = "wall_" + bc_name.rightJustified(3, '0');
          BoundaryCondition sym_bc(bc_name, "wall");
          GuiMainWindow::pointer()->addBC(bc, sym_bc);
        }
        set_bc();
        ++bc;
      }
    } while (!done);
  } else {
    if (0 <= mainWindow()->getPickedCell() && mainWindow()->getPickedCell() < GuiMainWindow::pointer()->getGrid()->GetNumberOfCells() ) {
      set_bc.setFeatureAngle(m_Ui.doubleSpinBoxFeatureAngle->value());
      set_bc.setBC(m_Ui.spinBoxBoundaryCode->value());

      set_bc.setProcessAll(m_ButtonGroup->button(1)->isChecked());
      set_bc.setSelectAllVisible(m_ButtonGroup->button(2)->isChecked());
      set_bc.setOnlyPickedCell(m_ButtonGroup->button(3)->isChecked());
      set_bc.setOnlyPickedCellAndNeighbours(m_ButtonGroup->button(4)->isChecked());

      cout << "GuiMainWindow::getPickedCell()=" << mainWindow()->getPickedCell() << endl;
      set_bc.setStart(mainWindow()->getPickedCell());
      set_bc();
    } else {
      EG_ERR_RETURN("Please select a cell first.");
    }
  }
}
