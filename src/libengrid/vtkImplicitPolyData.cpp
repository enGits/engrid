// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
/*=========================================================================

  Program:   Boolean
  Module:    $RCSfile: vtkImplicitPolyData.cpp,v $
  Language:  C++
  Date:      $Date: 2003/07/12 11:22:38 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2002 Denis Shamonin
  Section Computational Science
  University of Amsterdam
  Kruislaan 403, 1098 SJ Amsterdam
  the Netherlands

  E-mail        : dshamoni@science.uva.nl
  URL           : http://www.science.uva.nl/~dshamoni/
  
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================
ToDo:
-   Needs to be restructured in line with vtk principles of on demand
    execution. Eg BuildLinks and BuildLocator on first call to Evaluate*?.
    Handle input modification correctly etc.
-   Drop internal triangle filter and just check cell type as encountered?.
-   Get working with CellLocator. Currently crashes every time, although
    a large number of cells are located successfully. Too hard to debug,
    maybe try a tiny dataset.

    PLEASE SEND ME YOUR UPDATES, BUG FIXES! at dshamoni@science.uva.nl
=========================================================================*/
//#include <stdafx.h>
#include "vtkImplicitPolyData.h"
#include "vtkPolygon.h"

// Constructor
vtkImplicitPolyData::vtkImplicitPolyData()
{
  this->NoGradient[0] = 0.0;
  this->NoGradient[1] = 0.0;
  this->NoGradient[2] = 1.0;

  this->tri = NULL;
  this->input = NULL;
  this->locator = NULL;
  this->poly = vtkPolygon::New();
  this->cells = NULL;
  this->EvaluateBoundsSet=0;
  this->Tolerance = 10.0;
}

void vtkImplicitPolyData::SetInput(vtkPolyData *input)
{
  if( this->input != input ) 
    {
    vtkDebugMacro( <<" setting Input to " << (void *)input );

    // use a tringle filter on the polydata input
    // this is done to filter out lines and vertices to leave only
    // polygons which are required by this algorithm for cell normals
    if( this->tri == NULL )
      {
      this->tri = vtkTriangleFilter::New();
      this->tri->PassVertsOff();
      this->tri->PassLinesOff();
      }
    this->tri->SetInput( input );
    this->tri->Update();

    this->input = this->tri->GetOutput();
    this->input->BuildLinks();	// to enable calls to GetPointCells
    this->NoValue = this->input->GetLength();

    if( this->locator != NULL ) this->locator->Delete();
    locator = PointLocator::New();
    this->locator->SetDataSet( this->input );
//  if( this->EvaluateBoundsSet )     // experimental
//  this->locator->SetEvaluateBounds( this->EvaluateBounds );
    this->locator->BuildLocator();

    if( this->cells   != NULL ) this->cells->Delete();
    this->cells = vtkIdList::New();
    this->cells->SetNumberOfIds( this->input->GetNumberOfCells() );

    this->poly = vtkPolygon::New();
    this->Modified();
    }
}

// used to extend bounds of point locator
void vtkImplicitPolyData::SetEvaluateBounds( double eBounds[6] )
{
  int i;
  for( i=0; i<6; i++ ) 
    {
    this->EvaluateBounds[i] = eBounds[i];
    }
    
  this->EvaluateBoundsSet = 1;
}

unsigned long vtkImplicitPolyData::GetMTime()
{
  unsigned long mTime=this->vtkImplicitFunction::GetMTime();
  unsigned long inputMTime;

  if ( this->input != NULL )
    {
    this->input->Update ();
    inputMTime = this->input->GetMTime();
    mTime = ( inputMTime > mTime ? inputMTime : mTime );
    }

  return mTime;
}

// Destructor
vtkImplicitPolyData::~vtkImplicitPolyData()
{
  if( this->tri     != NULL ) this->tri->Delete();
  if( this->locator != NULL ) this->locator->Delete();
  if( this->poly    != NULL ) this->poly->Delete();
  if( this->cells   != NULL ) this->cells->Delete();
}

// Evaluate for point x[3].
// Method using combination of implicit planes
double vtkImplicitPolyData::EvaluateFunction(double x[3])
{
  // See if data set with polygons has been specified
  if( this->input == NULL || input->GetNumberOfCells() == 0 ) 
    {
    vtkErrorMacro(<<"No polygons to evaluate function!");
    return this->NoValue;
    }

  int cellNum, pid;
  double dist;
  vtkCell *cell;
  double dot, ret=-VTK_LARGE_FLOAT, cNormal[3], closestPoint[3];

  // get point id of closest point in data set according Tolerance
  pid = this->locator->FindClosestPointWithinRadius(Tolerance,x,dist);
  if(( pid != -1 ) && (dist < Tolerance))
    {
    this->input->GetPoint( pid, closestPoint );

    // get cells it belongs to
    this->input->GetPointCells( pid, cells );
    // for each cell
    for( cellNum=0; cellNum<cells->GetNumberOfIds(); cellNum++ )
      {
      cell = this->input->GetCell( cells->GetId( cellNum ) );

      // get cell normal
      poly->ComputeNormal( cell->GetPoints(), cNormal );

      dot = ( cNormal[0]*(x[0]-closestPoint[0]) +
              cNormal[1]*(x[1]-closestPoint[1]) +
              cNormal[2]*(x[2]-closestPoint[2]) );

      if( dot > ret ) ret = dot;
      }
    }
    
  if( ret == -VTK_LARGE_FLOAT ) ret = NoValue;
  return ret;
}

// Evaluate function gradient at point x[3].
void vtkImplicitPolyData::EvaluateGradient( double x[3], double n[3] )
{
  int i;
  // See if data set with polygons has been specified
  if( this->input == NULL || input->GetNumberOfCells() == 0 )
    {
    vtkErrorMacro(<<"No polygons to evaluate gradient!");
    for( i=0; i<3; i++ ) n[i] = this->NoGradient[i];
    return;
    }
    
  int cellNum, pid;
  double dist;
  vtkCell *cell;
  double dot, ret=-VTK_LARGE_FLOAT, cNormal[3], closestPoint[3];
  
  // get point id of closest point in data set according Tolerance
  pid = this->locator->FindClosestPointWithinRadius(Tolerance,x,dist);
  if(( pid != -1 ) && (dist < Tolerance))
    {
    this->input->GetPoint( pid, closestPoint );

    // get cells it belongs to
    this->input->GetPointCells( pid, cells );
    // for each cell
    for( cellNum=0; cellNum<cells->GetNumberOfIds(); cellNum++ ) 
       {
       cell = this->input->GetCell( cells->GetId( cellNum ) );

       // get cell normal
       poly->ComputeNormal( cell->GetPoints(), cNormal );

       dot = ( cNormal[0]*(x[0]-closestPoint[0]) +
               cNormal[1]*(x[1]-closestPoint[1]) +
               cNormal[2]*(x[2]-closestPoint[2]) );

       if( dot > ret ) 
         {
         for( i=0; i<3; i++ ) n[i] = cNormal[i];
         }
	 
	}
    }
    
  if( ret == -VTK_LARGE_FLOAT ) 
    {
    for( i=0; i<3; i++ ) n[i] = this->NoGradient[i];
    }
}

void vtkImplicitPolyData::PrintSelf(ostream& os, vtkIndent indent)
{
  vtkImplicitFunction::PrintSelf(os,indent);

  os << indent << "No polydata Value: " << this->NoValue << "\n";
  os << indent << "No polydata Gradient: (" << this->NoGradient[0] << ", "
     << this->NoGradient[1] << ", " << this->NoGradient[2] << ")\n";

  if ( this->input )
    {
    os << indent << "Input : " << this->input << "\n";
    }
  else
    {
    os << indent << "Input : (none)\n";
    }
  
  os << indent << "Tolerance: " << this->Tolerance << "\n";    
}
