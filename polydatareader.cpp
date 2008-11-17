#include "polydatareader.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"

#include <vtkXMLPolyDataReader.h>

PolyDataReader::PolyDataReader()
{
  setFormat("VTK poly data files(*.vtp)");
  setExtension(".vtp");
};

void PolyDataReader::operate()
{
  try {
    readInputFileName();
    if (isValid()) {
      EG_VTKSP(vtkXMLPolyDataReader,vtp);
      EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter,pd2ug);
      vtp->SetFileName(getFileName().toAscii().data());
      pd2ug->SetInput(vtp->GetOutput());
      pd2ug->Update();
      grid->DeepCopy(pd2ug->GetOutput());
      createBasicFields(grid, grid->GetNumberOfCells(), grid->GetNumberOfPoints(), false);
      UpdateNodeIndex(grid);
      UpdateCellIndex(grid);
      EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
      for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfPoints(); ++id_cell) {
        bc->SetValue(id_cell,99);
      };
    };
  } catch (Error err) {
    err.display();
  };
};
