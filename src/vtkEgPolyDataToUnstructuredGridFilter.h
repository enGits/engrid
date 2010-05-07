//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#ifndef __vtkEgPolyDataToUnstructuredGridFilter_h
#define __vtkEgPolyDataToUnstructuredGridFilter_h

class vtkEgPolyDataToUnstructuredGridFilter;

#include <vtkUnstructuredGridAlgorithm.h>
#include "engrid.h"

class vtkEgPolyDataToUnstructuredGridFilter : public vtkUnstructuredGridAlgorithm
{
  
public:
  
  static vtkEgPolyDataToUnstructuredGridFilter* New();
  
protected:
  
  vtkEgPolyDataToUnstructuredGridFilter();
  ~vtkEgPolyDataToUnstructuredGridFilter();
  
  virtual int FillInputPortInformation
    (
      int, 
      vtkInformation* info
    );
  
  virtual int FillOutputPortInformation
    (
      int, 
      vtkInformation* info
    );
  
  virtual int RequestData
    (
      vtkInformation*, 
      vtkInformationVector** inputVector, 
      vtkInformationVector*  outputVector
    );
  
private:
  
  vtkEgPolyDataToUnstructuredGridFilter (const vtkEgPolyDataToUnstructuredGridFilter&);
  void operator= (const vtkEgPolyDataToUnstructuredGridFilter&);
  
};

#endif
