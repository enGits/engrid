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

void VtkReader::createBoundaryFaces()
{
  QList<QVector<vtkIdType> > all_faces;
  MeshPartition part(m_Grid, true);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    for (int i = 0; i < part.c2cGSize(id_cell); ++i) {
      if (part.c2cGG(id_cell, i) < 0) {
        QVector<vtkIdType> face;
        getFaceOfCell(m_Grid, id_cell, i, face);
        all_faces << face;
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, m_Grid->GetNumberOfCells() + all_faces.size(), m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    copyCell(m_Grid, id_cell, grid);
  }
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  foreach (QVector<vtkIdType> face, all_faces) {
    EG_VTKSP(vtkIdList, pts);
    pts->SetNumberOfIds(face.size());
    for (int i = 0; i < face.size(); ++i) {
      pts->SetId(i, face[i]);
    }
    vtkIdType id_cell;
    if (face.size() == 3) {
      id_cell = grid->InsertNextCell(VTK_TRIANGLE, pts);
    } else if (face.size() == 4) {
      id_cell = grid->InsertNextCell(VTK_QUAD, pts);
    } else {
      id_cell = grid->InsertNextCell(VTK_POLYGON, pts);
    }
    bc->SetValue(id_cell, 1);
  }
  makeCopy(grid, m_Grid);
  GuiMainWindow::pointer()->setBC(1, BoundaryCondition("all_faces", "patch", 1));
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
      //m_Grid->DeepCopy(vtk->GetOutput());
      makeCopy(vtk->GetOutput(), m_Grid);
      createBoundaryFaces();
      createBasicFields(m_Grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
      UpdateNodeIndex(m_Grid);
      UpdateCellIndex(m_Grid);
    }
  } catch (Error err) {
    err.display();
  }
}
