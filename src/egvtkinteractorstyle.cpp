//
// C++ Implementation: egvtkinteractorstyle
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "egvtkinteractorstyle.h"

#include "vtkInteractorStyleUser.h"
#include "vtkMath.h"
#include "vtkCellPicker.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h"
#include "vtkCommand.h"

vtkCxxRevisionMacro(egvtkInteractorStyle, "$Revision: 1.37 $");
vtkStandardNewMacro(egvtkInteractorStyle);

//----------------------------------------------------------------------------
egvtkInteractorStyle::egvtkInteractorStyle() 
{
  this->MotionFactor   = 10.0;
}

//----------------------------------------------------------------------------
egvtkInteractorStyle::~egvtkInteractorStyle() 
{
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnMouseMove() 
{ 
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];
  
  switch (this->State) 
  {
  case VTKIS_ROTATE:
    this->FindPokedRenderer(x, y);
    this->Rotate();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;
    
  case VTKIS_PAN:
    this->FindPokedRenderer(x, y);
    this->Pan();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;
    
  case VTKIS_DOLLY:
    this->FindPokedRenderer(x, y);
    this->Dolly();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;
    
  case VTKIS_SPIN:
    this->FindPokedRenderer(x, y);
    this->Spin();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;
  }
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnLeftButtonDown()
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetShiftKey()) 
  {
    if (this->Interactor->GetControlKey()) 
    {
      this->StartDolly();
    }
    else 
    {
      this->StartPan();
    }
  } 
  else 
  {
    if (this->Interactor->GetControlKey())
    {
      this->StartSpin();
    }
    else 
    {
      this->StartRotate();
    }
  }
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnLeftButtonUp()
{
  switch (this->State) 
  {
  case VTKIS_DOLLY:
    this->EndDolly();
    break;
    
  case VTKIS_PAN:
    this->EndPan();
    break;
    
  case VTKIS_SPIN:
    this->EndSpin();
    break;
    
  case VTKIS_ROTATE:
    this->EndRotate();
    break;
  }
  
  if ( this->Interactor )
  {
    this->ReleaseFocus();
  }
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnMiddleButtonDown() 
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  this->GrabFocus(this->EventCallbackCommand);
  this->StartPan();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnMiddleButtonUp()
{
  switch (this->State) 
  {
  case VTKIS_PAN:
    this->EndPan();
    if ( this->Interactor )
    {
      this->ReleaseFocus();
    }
    break;
  }
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnRightButtonDown() 
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];
  cout<<"You clicked at ("<<X<<","<<Y<<")"<<endl;
  
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  this->GrabFocus(this->EventCallbackCommand);
  this->StartDolly();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnRightButtonUp()
{
  switch (this->State) 
  {
  case VTKIS_DOLLY:
    this->EndDolly();
    
    if ( this->Interactor )
    {
      this->ReleaseFocus();
    }
    break;
  }
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnMouseWheelForward() 
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  this->GrabFocus(this->EventCallbackCommand);
  this->StartDolly();
  double factor = this->MotionFactor * 0.2 * this->MouseWheelMotionFactor;
  this->Dolly(pow(1.1, factor));
  this->EndDolly();
  this->ReleaseFocus();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::OnMouseWheelBackward()
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  this->GrabFocus(this->EventCallbackCommand);
  this->StartDolly();
  double factor = this->MotionFactor * -0.2 * this->MouseWheelMotionFactor;
  this->Dolly(pow(1.1, factor));
  this->EndDolly();
  this->ReleaseFocus();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::Rotate()
{
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  
  int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  
  int *size = this->CurrentRenderer->GetRenderWindow()->GetSize();
  
  double delta_elevation = -20.0 / size[1];
  double delta_azimuth = -20.0 / size[0];
  
  double rxf = dx * delta_azimuth * this->MotionFactor;
  double ryf = dy * delta_elevation * this->MotionFactor;
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->Azimuth(rxf);
  camera->Elevation(ryf);
  camera->OrthogonalizeViewUp();
  
  if (this->AutoAdjustCameraClippingRange)
  {
    this->CurrentRenderer->ResetCameraClippingRange();
  }
  
  if (rwi->GetLightFollowCamera()) 
  {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  }
  
  rwi->Render();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::Spin()
{
  if ( this->CurrentRenderer == NULL )
  {
    return;
  }
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  
  double *center = this->CurrentRenderer->GetCenter();
  
  double newAngle = atan2( rwi->GetEventPosition()[1] - center[1], rwi->GetEventPosition()[0] - center[0] ) * vtkMath::RadiansToDegrees();
  
  double oldAngle = atan2( rwi->GetLastEventPosition()[1] - center[1], rwi->GetLastEventPosition()[0] - center[0] ) * vtkMath::RadiansToDegrees();
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->Roll( newAngle - oldAngle );
  camera->OrthogonalizeViewUp();
  
  rwi->Render();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::Pan()
{
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  
  double viewFocus[4], focalDepth, viewPoint[3];
  double newPickPoint[4], oldPickPoint[4], motionVector[3];
  
  // Calculate the focal depth since we'll be using it a lot
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->GetFocalPoint(viewFocus);
  this->ComputeWorldToDisplay(viewFocus[0], viewFocus[1], viewFocus[2], 
                              viewFocus);
  focalDepth = viewFocus[2];
  
  this->ComputeDisplayToWorld(rwi->GetEventPosition()[0], 
                              rwi->GetEventPosition()[1],
                              focalDepth, 
                              newPickPoint);
  
  // Has to recalc old mouse point since the viewport has moved,
  // so can't move it outside the loop
  
  this->ComputeDisplayToWorld(rwi->GetLastEventPosition()[0],
                              rwi->GetLastEventPosition()[1],
                              focalDepth, 
                              oldPickPoint);
  
  // Camera motion is reversed
  
  motionVector[0] = oldPickPoint[0] - newPickPoint[0];
  motionVector[1] = oldPickPoint[1] - newPickPoint[1];
  motionVector[2] = oldPickPoint[2] - newPickPoint[2];
  
  camera->GetFocalPoint(viewFocus);
  camera->GetPosition(viewPoint);
  camera->SetFocalPoint(motionVector[0] + viewFocus[0],
                        motionVector[1] + viewFocus[1],
                        motionVector[2] + viewFocus[2]);
  
  camera->SetPosition(motionVector[0] + viewPoint[0],
                      motionVector[1] + viewPoint[1],
                      motionVector[2] + viewPoint[2]);
  
  if (rwi->GetLightFollowCamera()) 
  {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  }
  
  rwi->Render();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::Dolly()
{
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  double *center = this->CurrentRenderer->GetCenter();
  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  double dyf = this->MotionFactor * dy / center[1];
  this->Dolly(pow(1.1, dyf));
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::Dolly(double factor)
{
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  
  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  if (camera->GetParallelProjection())
  {
    camera->SetParallelScale(camera->GetParallelScale() / factor);
  }
  else
  {
    camera->Dolly(factor);
    if (this->AutoAdjustCameraClippingRange)
    {
      this->CurrentRenderer->ResetCameraClippingRange();
    }
  }
  
  if (this->Interactor->GetLightFollowCamera()) 
  {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  }
  
  this->Interactor->Render();
}

//----------------------------------------------------------------------------
void egvtkInteractorStyle::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "MotionFactor: " << this->MotionFactor << "\n";
}

void egvtkInteractorStyle::OnChar()
{
  cout<<"OnChar "<<this->Interactor->GetKeyCode()<<endl;
  this->EventCallbackCommand->SetAbortFlag(1);
  
  vtkRenderWindowInteractor *rwi = this->Interactor;
  
  switch (rwi->GetKeyCode()) 
  {
    case 'n' :
      cout<<"pick node by mouse"<<endl;
      break;
    case 'N' :
      cout<<"pick node by ID"<<endl;
      break;
    
    case 'c' :
      cout<<"pick cell by mouse"<<endl;
      break;
    case 'C' :
      cout<<"pick cell by ID"<<endl;
      break;
    
    case 'b' :
      cout<<"box select"<<endl;
    
      break;
  }
  
   // otherwise pass the OnChar to the vtkInteractorStyle.
  if (this->HasObserver(vtkCommand::CharEvent)) 
    {
    this->ShiftKey = this->Interactor->GetShiftKey();
    this->ControlKey = this->Interactor->GetControlKey();
    this->KeyCode = this->Interactor->GetKeyCode();  
    
    this->InvokeEvent(vtkCommand::CharEvent,NULL);
    }
  else
    {
    this->vtkInteractorStyle::OnChar();
    }
  
}
