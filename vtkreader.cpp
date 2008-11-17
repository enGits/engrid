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
