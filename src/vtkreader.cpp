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
#include "vtkreader.h"

#include <vtkUnstructuredGridReader.h>

VtkReader::VtkReader()
{
  setFormat("legacy VTK files(*.vtk)");
  setExtension(".vtk");
};

void VtkReader::operate()
{
  try {
    readInputFileName();
    if (isValid()) {
      EG_VTKSP(vtkUnstructuredGridReader,vtk);
      vtk->SetFileName(getFileName().toAscii().data());
      vtk->Update();
      grid->DeepCopy(vtk->GetOutput());
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints(), false);
      UpdateNodeIndex(grid);
      UpdateCellIndex(grid);
    };
  } catch (Error err) {
    err.display();
  };
};
