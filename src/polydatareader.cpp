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
#include "polydatareader.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"

#include <vtkXMLPolyDataReader.h>

PolyDataReader::PolyDataReader()
{
  setFormat("VTK poly data files(*.vtp)");
  setExtension(".vtp");
}

void PolyDataReader::operate()
{
  try {
    readInputFileName();
    if (isValid()) {
      EG_VTKSP(vtkXMLPolyDataReader,vtp);
      EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter,pd2ug);
      vtp->SetFileName(qPrintable(getFileName()));
      pd2ug->SetInput(vtp->GetOutput());
      pd2ug->Update();
      grid->DeepCopy(pd2ug->GetOutput());
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints());
      UpdateNodeIndex(grid);
      UpdateCellIndex(grid);
      EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
      for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfPoints(); ++id_cell) {
        bc->SetValue(id_cell,99);
      }
    }
  } catch (Error err) {
    err.display();
  }
}
