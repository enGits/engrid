//
// C++ Implementation: boxselect
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
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
/*      vtkTransform *t = vtkTransform::New();
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      widget->GetTransform(t);
      widget->GetProp3D()->SetUserTransform(t);
      t->Delete();*/
  }
};

BoxSelect::BoxSelect()
 : Operation()
{
}


BoxSelect::~BoxSelect()
{
}


void BoxSelect::operate()
{
  cout<<"BoxSelect"<<endl;
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
  
}
