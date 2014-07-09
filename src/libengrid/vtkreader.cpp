// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "vtkreader.h"

#include <vtkUnstructuredGridReader.h>

#include <QFileInfo>
#include "guimainwindow.h"

VtkReader::VtkReader()
{
  setFormat("legacy VTK files(*.vtk)");
  setExtension(".vtk");
}

void VtkReader::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readInputFileName(file_info.completeBaseName() + ".vtk");
    if (isValid()) {
      EG_VTKSP(vtkUnstructuredGridReader,vtk);
      vtk->SetFileName(qPrintable(getFileName()));
      vtk->Update();
      m_Grid->DeepCopy(vtk->GetOutput());
      createBasicFields(m_Grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
      UpdateNodeIndex(m_Grid);
      UpdateCellIndex(m_Grid);
    }
  } catch (Error err) {
    err.display();
  }
}
