//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "vtkEgPolyDataToUnstructuredGridFilter.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(vtkEgPolyDataToUnstructuredGridFilter);

vtkEgPolyDataToUnstructuredGridFilter::vtkEgPolyDataToUnstructuredGridFilter()
  : vtkUnstructuredGridAlgorithm()
{
  this->SetNumberOfInputPorts(1);
};

vtkEgPolyDataToUnstructuredGridFilter::~vtkEgPolyDataToUnstructuredGridFilter()
{
};

int vtkEgPolyDataToUnstructuredGridFilter::FillInputPortInformation
(
  int, 
  vtkInformation* info
)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkPolyData");
  return 1;
};

int vtkEgPolyDataToUnstructuredGridFilter::FillOutputPortInformation
(
  int, 
  vtkInformation* info
)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkUnstructuredGrid");
  return 1;
};

int vtkEgPolyDataToUnstructuredGridFilter::RequestData
(
  vtkInformation*, 
  vtkInformationVector** inputVector, 
  vtkInformationVector*  outputVector
)
{
  try {
    
    // Get input and output data
    vtkPolyData*         input  = vtkPolyData::GetData(inputVector[0]);
    vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData(outputVector);
    
    // the coordinates of the nodes
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(input->GetPoints()->GetNumberOfPoints());
    for (vtkIdType vertId = 0; vertId < input->GetPoints()->GetNumberOfPoints(); ++vertId) {
      double x[3];
      input->GetPoints()->GetPoint(vertId,x);
      points->SetPoint(vertId,x);
    };
    output->SetPoints(points);
    
    // the cells
    input->BuildCells();
    output->Allocate(input->GetNumberOfPolys(),input->GetNumberOfPolys());
    for (vtkIdType cellId = 0; cellId < input->GetNumberOfPolys(); ++cellId) {
      vtkIdType *pts;
      vtkIdType npts;
      input->GetCellPoints(cellId,npts,pts);
      if        (npts == 3) {
        output->InsertNextCell(VTK_TRIANGLE,3,pts);
      } else if (npts == 4) {
        output->InsertNextCell(VTK_QUAD,4,pts);
      };
    };
  
  } catch (Error err) {
    err.display();
  };
  
  return 1;
  
};

