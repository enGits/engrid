#include "guipick.h"
#include "guimainwindow.h"

void GuiPick::before()
{
  ui.spinBox_Point->setMaximum(grid->GetNumberOfPoints()-1);
  ui.spinBox_Cell->setMaximum(grid->GetNumberOfCells()-1);
}

void GuiPick::operate()
{
  cout<<"GuiPick called"<<endl;
  if(ui.radioButton_Point->isChecked())
  {
    vtkIdType nodeId=ui.spinBox_Point->value();
    cout<<"Pick point "<<nodeId<<endl;
    GuiMainWindow::pointer()->setPickMode(false,false);
    GuiMainWindow::pointer()->pickPoint(nodeId);
  }
  else
  {
    vtkIdType cellId=ui.spinBox_Cell->value();
    cout<<"Pick cell "<<cellId<<endl;
    GuiMainWindow::pointer()->setPickMode(false,true);
    GuiMainWindow::pointer()->pickCell(cellId);
  }
  updateActors();
};
