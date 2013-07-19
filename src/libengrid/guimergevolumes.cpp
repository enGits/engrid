// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "guimergevolumes.h"
#include "guimainwindow.h"

GuiMergeVolumes::GuiMergeVolumes()
{
}

void GuiMergeVolumes::before()
{
  populateVolumes(m_Ui.listWidgetVC1);
  populateVolumes(m_Ui.listWidgetVC2);
}

void GuiMergeVolumes::operate()
{
  QString vol_name1 = getSelectedVolume(m_Ui.listWidgetVC1);
  QString vol_name2 = getSelectedVolume(m_Ui.listWidgetVC2);
  VolumeDefinition V1 = GuiMainWindow::pointer()->getVol(vol_name1);
  VolumeDefinition V2 = GuiMainWindow::pointer()->getVol(vol_name2);
  QSet<int> all_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();

  // identify boundary patches to be deleted
  QSet<int> del_bcs;
  foreach (int bc, all_bcs) {
    int sign1 = V1.getSign(bc);
    int sign2 = V2.getSign(bc);
    if (sign1 != 0 && sign2 != 0) {
      if (sign1*sign2 > 0) {
        EG_ERR_RETURN("volume definition not consistent (green/yellow)");
      }
      del_bcs.insert(bc);
    }
  }

  // count cells of new mesh
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  int num_new_cells = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      ++num_new_cells;
    } else {
      if (!del_bcs.contains(cell_code->GetValue(id_cell))) {
        ++num_new_cells;
      }
    }
  }

  // copy nodes and cells
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, num_new_cells, m_Grid->GetNumberOfPoints());
  vtkIdType id_new_node = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_new_node, x.data());
    copyNodeData(m_Grid, id_node, new_grid, id_new_node);
    ++id_new_node;
  }
  vtkIdType id_new_cell;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    //vtkIdType N_pts, *pts;
    //m_Grid->GetCellPoints(id_cell, N_pts, pts);
    bool insert_cell = true;
    if (isSurface(id_cell, m_Grid)) {
      if (del_bcs.contains(cell_code->GetValue(id_cell))) {
        insert_cell = false;
      }
    }
    if (insert_cell) {
      id_new_cell = copyCell(m_Grid, id_cell, new_grid);
      //new_grid->InsertNextCell(m_Grid->GetCellType(id_cell), N_pts, pts);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    }
  }

  // update volume definitions
  QList<VolumeDefinition> old_vols = GuiMainWindow::pointer()->getAllVols();
  QList<VolumeDefinition> new_vols;
  foreach (VolumeDefinition V, old_vols) {
    if (V.getName() != vol_name1 && V.getName() != vol_name2) {
      new_vols.append(V);
    }
  }
  VolumeDefinition V(m_Ui.lineEditNewVol->text(), V1.getVC());
  foreach (int bc, all_bcs) {
    if (!del_bcs.contains(bc)) {
      if (V1.getSign(bc) != 0) {
        V.addBC(bc, V1.getSign(bc));
      } else if (V2.getSign(bc) != 0) {
        V.addBC(bc, V2.getSign(bc));
      } else {
        V.addBC(bc, 0);
      }
    }
  }
  new_vols.append(V);
  GuiMainWindow::pointer()->setAllVols(new_vols);
  makeCopy(new_grid, m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(false);
}
















