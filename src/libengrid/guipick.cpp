// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "guipick.h"
#include "guimainwindow.h"

void GuiPick::before()
{
  m_Ui.spinBox_Point->setMaximum(m_Grid->GetNumberOfPoints()-1);
  m_Ui.spinBox_Cell->setMaximum(m_Grid->GetNumberOfCells()-1);
}

void GuiPick::operate()
{
  cout<<"GuiPick called"<<endl;
  if(m_Ui.radioButton_Point->isChecked())
  {
    vtkIdType nodeId=m_Ui.spinBox_Point->value();
    cout<<"Pick point "<<nodeId<<endl;
    GuiMainWindow::pointer()->setPickMode(false,false);
    GuiMainWindow::pointer()->pickPoint(nodeId);
  }
  else
  {
    vtkIdType cellId=m_Ui.spinBox_Cell->value();
    cout<<"Pick cell "<<cellId<<endl;
    GuiMainWindow::pointer()->setPickMode(false,true);
    GuiMainWindow::pointer()->pickCell(cellId);
  }
  updateActors();
};
