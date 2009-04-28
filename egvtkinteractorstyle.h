#ifndef EGVTKINTERACTORSTYLE_H
#define EGVTKINTERACTORSTYLE_H

#include <vtkInteractorStyle.h>

#include "vtkCamera.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"

#include <iostream>
using namespace std;

class VTK_RENDERING_EXPORT egvtkInteractorStyle : public vtkInteractorStyle
{
public:
  static egvtkInteractorStyle *New();
  vtkTypeRevisionMacro(egvtkInteractorStyle,vtkInteractorStyle);
  void PrintSelf(ostream& os, vtkIndent indent);
  
//   vtkGetMacro(ShiftKey,int);
//   vtkGetMacro(CtrlKey,int);
  
  // Description:
  // Event bindings controlling the effects of pressing mouse buttons
  // or moving the mouse.
  virtual void OnMouseMove();
  virtual void OnLeftButtonDown();
  virtual void OnLeftButtonUp();
  virtual void OnMiddleButtonDown();
  virtual void OnMiddleButtonUp();
  virtual void OnRightButtonDown();
  virtual void OnRightButtonUp();
  virtual void OnMouseWheelForward();
  virtual void OnMouseWheelBackward();
  virtual void OnChar();
  
/*  virtual void 	OnKeyDown (){
    cout<<"OnKeyDown "<<this->Interactor->GetKeyCode()<<endl;
    this->EventCallbackCommand->SetAbortFlag(1);
  };
  virtual void 	OnKeyUp (){
    cout<<"OnKeyUp "<<this->Interactor->GetKeyCode()<<endl;
    this->EventCallbackCommand->SetAbortFlag(1);
  };
  virtual void 	OnKeyPress (){
    cout<<"OnKeyPress "<<this->Interactor->GetKeyCode()<<endl;
    this->EventCallbackCommand->SetAbortFlag(1);
  };
  virtual void 	OnKeyRelease (){
    cout<<"OnKeyRelease "<<this->Interactor->GetKeyCode()<<endl;
    this->EventCallbackCommand->SetAbortFlag(1);
  };*/
  
  // These methods for the different interactions in different modes
  // are overridden in subclasses to perform the correct motion. Since
  // they are called by OnTimer, they do not have mouse coord parameters
  // (use interactor's GetEventPosition and GetLastEventPosition)
  virtual void Rotate();
  virtual void Spin();
  virtual void Pan();
  virtual void Dolly();
  
  // Description:
  // Set the apparent sensitivity of the interactor style to mouse motion.
  vtkSetMacro(MotionFactor,double);
  vtkGetMacro(MotionFactor,double);
  
protected:
  egvtkInteractorStyle();
  ~egvtkInteractorStyle();
  
  double MotionFactor;
  
  virtual void Dolly(double factor);
  
private:
  egvtkInteractorStyle(const egvtkInteractorStyle&);  // Not implemented.
  void operator=(const egvtkInteractorStyle&);  // Not implemented.
};

#endif
