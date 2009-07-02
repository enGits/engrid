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
  EG_VTKDCC(vtkIntArray,cell_code,m_Input,"cell_code");
  
  // copy the node coordinates
  allocateGrid(m_Output, m_Input->GetNumberOfCells(), m_Input->GetNumberOfPoints());
  for (vtkIdType vertId = 0; vertId < m_Input->GetNumberOfPoints(); ++vertId) {
    double x[3];
    m_Input->GetPoints()->GetPoint(vertId,x);
    m_Output->GetPoints()->SetPoint(vertId,x);
  }
  
  // copy the cells and the cell/node data
  for (vtkIdType cellId = 0; cellId < m_Input->GetNumberOfCells(); ++cellId) {
    vtkIdType *pts;
    vtkIdType npts;
    bool add = false;
    if (!cell_code) {
      add = false;
    }
    else {
      if (m_BoundaryCodes.contains(cell_code->GetValue(cellId))) {
        if (m_Input->GetCellType(cellId) == VTK_TRIANGLE) add = true;
        if (m_Input->GetCellType(cellId) == VTK_QUAD)     add = true;
      }
    }
    
    if (add) {
      m_Input->GetCellPoints(cellId,npts,pts);
      vtkIdType newCell = m_Output->InsertNextCell(m_Input->GetCellType(cellId),npts,pts);
      copyCellData(m_Input, cellId, m_Output, newCell);
      for(int i = 0; i < npts; i++) {
        copyNodeData(m_Input, pts[i], m_Output, pts[i]);
      }
    }
  }
  
  m_Output->Squeeze();
}

