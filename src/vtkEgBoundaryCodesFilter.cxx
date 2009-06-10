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
#include "vtkEgBoundaryCodesFilter.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>

vtkStandardNewMacro(vtkEgBoundaryCodesFilter)

void vtkEgBoundaryCodesFilter::ExecuteEg()
{
  // the coordinates of the nodes
  allocateGrid(output, input->GetNumberOfCells(), input->GetNumberOfPoints());
  for (vtkIdType vertId = 0; vertId < input->GetNumberOfPoints(); ++vertId) {
    double x[3];
    input->GetPoints()->GetPoint(vertId,x);
    output->GetPoints()->SetPoint(vertId,x);
    copyNodeData(input, vertId, output, vertId);
  };
  
  EG_VTKDCC(vtkIntArray,cell_code,input,"cell_code");
  
  // the cells
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    vtkIdType *pts;
    vtkIdType npts;
    bool add = false;
    if (!cell_code || !BoundaryCodes) { 
      add = false;
    } else if (BoundaryCodes->contains(cell_code->GetValue(cellId))) {
      if (input->GetCellType(cellId) == VTK_TRIANGLE) add = true;
      if (input->GetCellType(cellId) == VTK_QUAD)     add = true;
    };
    if (add) {
      input->GetCellPoints(cellId,npts,pts);
      vtkIdType newCell = output->InsertNextCell(input->GetCellType(cellId),npts,pts);
      copyCellData(input, cellId, output, newCell);
    };
  };
  
  output->Squeeze();
};

