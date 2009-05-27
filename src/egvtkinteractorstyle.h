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
  
    // Event information
  int   AltKey;
  int   ControlKey;
  int   ShiftKey;
  char  KeyCode;
  int   RepeatCount;
  char* KeySym; 
  int   EventPosition[2];
  int   LastEventPosition[2];
  int   EventSize[2];
  int   Size[2];
  int   TimerEventId;
  int   TimerEventType;
  int   TimerEventDuration;
  int   TimerEventPlatformId;
  
  vtkSetMacro(AltKey, int);
  vtkGetMacro(AltKey, int);
  vtkSetMacro(ControlKey, int);
  vtkGetMacro(ControlKey, int);
  vtkSetMacro(ShiftKey, int);
  vtkGetMacro(ShiftKey, int);
  vtkSetMacro(KeyCode, char);
  vtkGetMacro(KeyCode, char);
  vtkSetMacro(RepeatCount, int);
  vtkGetMacro(RepeatCount, int);
  vtkSetStringMacro(KeySym);
  vtkGetStringMacro(KeySym);
  
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
