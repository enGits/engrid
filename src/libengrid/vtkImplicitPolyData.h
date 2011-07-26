// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
  Module:    $RCSfile: vtkImplicitPolyData.h,v $
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

=========================================================================*/
// .NAME vtkImplicitPolyData - treat polygons in input as implicit planes
// .SECTION Description
// vtkImplicitPolyData provides the basis for using arbritrary polygonal data
// as an implicit surface, points are evaluated against nearest polygons which
// are handled as implicit planes.
// vtkImplicitPolyData is a concrete implementation of the abstract class
// vtkImplicitFunction.
// An internal instance of vtkTriangleFilter is used to filter vertices and
// lines out of the input PolyData, and a MYvtkPointLocator is used to find the
// nearest triangle to a candidate point.

// PLEASE SEND ME YOUR UPDATES, BUG FIXES! at dshamoni@science.uva.nl

#ifndef __vtkImplicitPolyData_h
#define __vtkImplicitPolyData_h

#include <math.h>
#include "vtkPolyData.h"
#include "vtkImplicitFunction.h"
#include "vtkTriangleFilter.h"
#include "vtkIdList.h"

#define PointLocator vtkPointLocator
#include "vtkPointLocator.h"

class vtkImplicitPolyData : public vtkImplicitFunction
{
public:
  static vtkImplicitPolyData *New() {return new vtkImplicitPolyData;};
  const char *GetClassName() {return "vtkImplicitPolyData";};
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkImplicitPolyData();
  ~vtkImplicitPolyData();

  // Description:
  // Return the MTime also considering the Input dependency.
  unsigned long GetMTime();

  void SetEvaluateBounds( double eBounds[6] );
  
  // Description
  // Evaluate plane equation of nearest triangle to point x[3].
  double EvaluateFunction(double x[3]);

  // Description
  // Evaluate function gradient of nearest triangle to point x[3].
  void EvaluateGradient(double x[3], double g[3]);

  // Description:
  // Set the input polydata used for the implicit function evaluation.
  // Passes input through an internal instance of vtkTriangleFilter to remove
  // vertices and lines, leaving only triangular polygons for evaluation as
  // implicit planes
  void SetInput(vtkPolyData *input);

  // Description:
  // Set / get the function value to use if no input polydata specified.
  vtkSetMacro(NoValue,double);
  vtkGetMacro(NoValue,double);

  // Description:
  // Set / get the function gradient to use if no input polydata specified.
  vtkSetVector3Macro(NoGradient,double);
  vtkGetVector3Macro(NoGradient,double);

  vtkSetMacro(Tolerance,double);
  vtkGetMacro(Tolerance,double);

protected:
  double NoValue;
  double NoGradient[3];

  vtkTriangleFilter *tri;
  vtkPolyData *input;
  vtkPointLocator *locator;
  vtkPolygon *poly;
  vtkIdList *cells;

  double EvaluateBounds[6];
  int   EvaluateBoundsSet;
  double Tolerance;
  
};
#endif


