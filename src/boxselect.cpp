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
#include "boxselect.h"

#include "guimainwindow.h"

#include <vtkBoxWidget.h>
#include <vtkCommand.h>

class vtkMyCallback : public vtkCommand
{
public:
  static vtkMyCallback *New() 
  { return new vtkMyCallback; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {
    caller = caller; ///@@@ get rid of warning
/*      vtkTransform *t = vtkTransform::New();
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      widget->GetTransform(t);
      widget->GetProp3D()->SetUserTransform(t);
      t->Delete();*/
  }
};

BoxSelect::BoxSelect()
{
  EG_TYPENAME;
}


BoxSelect::~BoxSelect()
{
}


void BoxSelect::operate()
{
  cout<<"BoxSelect"<<endl;
  
///@@@ TODO: Finish BoxSelect
/*  
  boxWidget = vtkBoxWidget::New();
  boxWidget->SetInteractor(GuiMainWindow::pointer()->getInteractor());
//   this->Interactor;
  boxWidget->SetPlaceFactor(1.25);
  
  //
  // Place the interactor initially. The input to a 3D widget is used to 
  // initially position and scale the widget. The EndInteractionEvent is
  // observed which invokes the SelectPolygons callback.
  //
  GuiMainWindow::pointer()->getRenderer();
  GuiMainWindow::pointer()->getInteractor();
//   boxWidget->SetProp3D(surface_actor);
  boxWidget->PlaceWidget();
  vtkMyCallback *callback = vtkMyCallback::New();
  boxWidget->AddObserver(vtkCommand::InteractionEvent, callback);
  boxWidget->On();
  */
}
