// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#include "guicreatevolumemesh.h"
#include "createvolumemesh.h"
#include "guimainwindow.h"

GuiCreateVolumeMesh::GuiCreateVolumeMesh()
{
}

void GuiCreateVolumeMesh::before()
{
  GuiMainWindow::pointer()->createDefaultVol();
  populateVolumes(m_Ui.listWidget);
}

void GuiCreateVolumeMesh::operate()
{
  QString volume_name = getSelectedVolume(m_Ui.listWidget);
  VolumeDefinition V = mainWindow()->getVol(volume_name);
  CreateVolumeMesh mesh_volume;
  EG_VTKSP(vtkUnstructuredGrid, part_grid);
  EG_VTKSP(vtkUnstructuredGrid, rest_grid);
  {
    MeshPartition part(volume_name);
    part.setVolumeOrientation();
    part.extractToVtkGrid(part_grid);
    MeshPartition rest(m_Grid);
    rest.setRemainder(part);
    rest.extractToVtkGrid(rest_grid);
  }
  mesh_volume.setGrid(part_grid);
  mesh_volume();
  EG_VTKDCC(vtkIntArray, cell_code, part_grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < part_grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, part_grid)) {
      cell_code->SetValue(id_cell, V.getVC());
    }
  }
  makeCopy(part_grid, m_Grid);
  {
    MeshPartition part(m_Grid, true);
    MeshPartition rest(rest_grid, true);
    part.addPartition(rest);
  }
  resetOrientation(m_Grid);
  createIndices(m_Grid);
}


