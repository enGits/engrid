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
#ifndef tvtkoperation_H
#define tvtkoperation_H


template <class T>
class TVtkOperation;

#include "operation.h"

template <class T>
class TVtkOperation : public Operation
{
  
protected: // attributes
  
  T *vtk_filter;
  
public: // methods
  
  TVtkOperation();
  virtual ~TVtkOperation();
  T* getFilter() { return vtk_filter; };
  
protected: // methods
  
  virtual void operate();
  
};


template <class T>
TVtkOperation<T>::TVtkOperation()
{
  vtk_filter = T::New();
};

template <class T>
TVtkOperation<T>::~TVtkOperation()
{
  vtk_filter->Delete();
};

template <class T>
void TVtkOperation<T>::operate()
{
  vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
  ug->DeepCopy(grid);
  vtk_filter->SetInput(ug);
  vtk_filter->Update();
  grid->DeepCopy(vtk_filter->GetOutput());
};


#endif
