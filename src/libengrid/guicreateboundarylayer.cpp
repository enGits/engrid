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
#include "guicreateboundarylayer.h"
#include "guimainwindow.h"
#include "createboundarylayershell.h"
#include "deletecells.h"
#include "meshpartition.h"

GuiCreateBoundaryLayer::GuiCreateBoundaryLayer()
{
}

void GuiCreateBoundaryLayer::before()
{
}

void GuiCreateBoundaryLayer::operate()
{
  if (!GuiMainWindow::pointer()->checkCadInterfaces()) {
    GuiMainWindow::pointer()->storeCadInterfaces();
  }

  resetOrientation(m_Grid);

  CreateBoundaryLayerShell blayer;
  blayer.setGrid(m_Grid);
  blayer.setAllCells();
  blayer.setReverseBoundaryFaces(m_Ui.checkBoxMarchInA->isChecked());
  blayer.setMarchingDirectionA(false);
  blayer();
  if (!blayer.success()) {
    EG_ERR_RETURN("unable to create boundary layer shell");
  }

  vtkUnstructuredGrid *prismatic_grid = blayer.getPrismaticGrid();
  MeshPartition prismatic_part(prismatic_grid, true);

  QVector<vtkIdType> shell_cells;
  getSurfaceCells(blayer.getBoundaryLayerCodes(), shell_cells, m_Grid);
  if (shell_cells.size() == 0) {
    EG_ERR_RETURN("unable to extract outer boundary layer faces");
  }
  MeshPartition shell_part(m_Grid);
  shell_part.setCells(shell_cells);
  EG_VTKSP(vtkUnstructuredGrid, shell_grid);
  shell_part.extractToVtkGrid(shell_grid);
  if (blayer.usePolyhedralBoundaryLayer() && blayer.triangulateFarfield()) {
    blayer.triangulateFarfieldSurface(shell_grid);
  }

  auto new_bc = GuiMainWindow::pointer()->getBC(BoundaryCondition("farfield", "patch"));
  EG_VTKDCC(vtkIntArray, cell_code, shell_grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, shell_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, shell_grid, "cell_curdir");
  for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
    cell_code->SetValue(id_cell, new_bc.getCode());
    reorientateFace(shell_grid, id_cell);
    cell_orgdir->SetValue(id_cell, cell_curdir->GetValue(id_cell));
  }

  DeleteCells delete_cells;
  MeshPartition all(m_Grid, true);
  delete_cells.setGrid(m_Grid);
  delete_cells.setCellsToDelete(all.getCells());
  delete_cells();

  m_Part.setAllCells();
  m_Part.addPartition(prismatic_part);
  MeshPartition new_shell_part(shell_grid, true);
  m_Part.addPartition(new_shell_part);

  resetOrientation(m_Grid);
  createIndices(m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
  GuiMainWindow::pointer()->updateActors();
}
