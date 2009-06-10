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
#include "vtkEgGridFilter.h"
#include "vtkInformation.h"

vtkEgGridFilter::vtkEgGridFilter()
{
  BoundaryCodes = NULL;
  input = NULL;
  output = NULL;
}

int vtkEgGridFilter::FillInputPortInformation
(
  int, 
  vtkInformation* info
)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkUnstructuredGrid");
  return 1;
}

int vtkEgGridFilter::FillOutputPortInformation
(
  int, 
  vtkInformation* info
)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkUnstructuredGrid");
  return 1;
}

int vtkEgGridFilter::RequestData
(
  vtkInformation*, 
  vtkInformationVector** inputVector, 
  vtkInformationVector*  outputVector
)
{
  // Get input and output data
  input  = vtkUnstructuredGrid::GetData(inputVector[0]);
  output = vtkUnstructuredGrid::GetData(outputVector);
  
  // check if we have a proper inputVector
  if (!input)                                      { vtkDebugMacro("No input!"); return 1; }
  if (!input->GetPoints())                         { vtkDebugMacro("No input!"); return 1; }
  if (input->GetPoints()->GetNumberOfPoints() < 1) { vtkDebugMacro("No input!"); return 1; }
  
  // call the enGrid filter method
  try {
    ExecuteEg();
  } catch (Error err) {
    err.display();
    output->DeepCopy(input);
  }
  
  return 1;
}

void vtkEgGridFilter::SetBoundaryCodes(QSet<int> *bc)
{
  BoundaryCodes = bc;
  Modified();
}

void vtkEgGridFilter::ExtractBoundary
(
  QVector<vtkIdType>  &cells, 
  QVector<vtkIdType>  &nodes, 
  QSet<int>           &bc,
  vtkUnstructuredGrid *grid
)
{
  cells.clear();
  nodes.clear();
  QSet<vtkIdType> ex_nodes, ex_cells;
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
    if ((grid->GetCellType(i) == VTK_TRIANGLE) || (grid->GetCellType(i) == VTK_QUAD)){
      if (bc.contains(cell_code->GetValue(i))) {
        ex_cells.insert(i);
        vtkIdType *pts;
        vtkIdType npts;
        input->GetCellPoints(i,npts,pts);
        for (int j = 0; j < npts; ++j) {
          ex_nodes.insert(pts[j]);
        }
      }
    }
  }
  cells.resize(ex_cells.size());
  nodes.resize(ex_nodes.size());
  {
    int j = 0;
    vtkIdType i;
    foreach(i,ex_cells) {
      cells[j] = i;
      ++j;
    }
  }
  {
    int j = 0;
    vtkIdType i;
    foreach(i,ex_nodes) {
      nodes[j] = i;
      ++j;
    }
  }
}

/*
void vtkEgGridFilter::CreateBasicFields(vtkIdType Ncells, vtkIdType Nnodes)
{
  Nnodes = Nnodes;
  if (!output->GetCellData()->GetArray("cell_code")) {
    EG_VTKSP(vtkIntArray, cell_code);
    cell_code->SetName("cell_code");
    cell_code->SetNumberOfValues(Ncells);
    output->GetCellData()->AddArray(cell_code);
  } else {
    EG_VTKDCC(vtkIntArray, cell_code, output, "cell_code");
    cell_code->SetNumberOfValues(Ncells);
  }
  if (!output->GetCellData()->GetArray("cell_index")) {
    EG_VTKSP(vtkLongArray_t, cell_index);
    cell_index->SetName("cell_index");
    cell_index->SetNumberOfValues(Ncells);
    output->GetCellData()->AddArray(cell_index);
  } else {
    EG_VTKDCC(vtkLongArray_t, cell_index, output, "cell_index");
    cell_index->SetNumberOfValues(Ncells);
  }
}

void vtkEgGridFilter::CopyCellData(vtkIdType oldId, vtkIdType newId)
{
  EG_VTKDCC(vtkIntArray, cell_code1, input,  "cell_code");
  EG_VTKDCC(vtkIntArray, cell_code2, output, "cell_code");
  cell_code2->SetValue(newId, cell_code1->GetValue(oldId));
  EG_VTKDCC(vtkLongArray_t, cell_index1, input,  "cell_index");
  EG_VTKDCC(vtkLongArray_t, cell_index2, output, "cell_index");
  cell_index2->SetValue(newId, cell_index1->GetValue(oldId));
  //qDebug() << oldId << newId;
}

void vtkEgGridFilter::AllocateOutput(vtkIdType Ncells, vtkIdType Nnodes)
{
  EG_VTKSP(vtkPoints,points);
  points->SetNumberOfPoints(Nnodes);
  output->SetPoints(points);
  output->Allocate(Ncells,max(1,Ncells/10));
  CreateBasicFields(Ncells, Nnodes);
}
*/























