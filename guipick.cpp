#include "guipick.h"
#include "guimainwindow.h"

void GuiPick::before()
{
  ui.spinBox_Point->setMaximum(grid->GetNumberOfPoints());
  ui.spinBox_Cell->setMaximum(grid->GetNumberOfCells());
}

void GuiPick::operate()
{
  cout<<"GuiPick called"<<endl;
  if(ui.radioButton_Point->isChecked())
  {
    vtkIdType nodeId=ui.spinBox_Point->value();
    cout<<"Pick point "<<nodeId<<endl;
    GuiMainWindow::pointer()->pickPoint(nodeId);
  }
  else
  {
    vtkIdType cellId=ui.spinBox_Cell->value();
    cout<<"Pick cell "<<cellId<<endl;
    GuiMainWindow::pointer()->pickCell(cellId);
  }
  updateActors();
};
