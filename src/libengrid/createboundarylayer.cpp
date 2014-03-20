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

#include "createboundarylayer.h"

#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"

CreateBoundaryLayer::CreateBoundaryLayer()
{
  m_RestGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
}

void CreateBoundaryLayer::prepare()
{
  /*
  // set m_Grid to selected volume
  m_VolDef = GuiMainWindow::pointer()->getVol(m_VolumeName);
  foreach (int bc, m_BoundaryCodes) {
    if (m_VolDef.getSign(bc) == 0) {
      QString msg;
      msg.setNum(bc);
      msg = "Boundary code " + msg + " is not part of the volume '" + m_VolumeName +"'.";
      EG_ERR_RETURN(msg);
    }
  }

  EG_VTKSP(vtkUnstructuredGrid, rest_grid);
  {
    EG_VTKSP(vtkUnstructuredGrid, vol_grid);
    MeshPartition volume(m_VolumeName);
    MeshPartition rest(m_Grid);
    rest.setRemainder(volume);
    volume.setVolumeOrientation();
    volume.extractToVtkGrid(vol_grid);
    rest.extractToVtkGrid(m_RestGrid);
    makeCopy(vol_grid, m_Grid);
  }
  */

  readSettings();

  setAllCells();

  getSurfaceCells(m_BoundaryCodes, layer_cells, m_Grid);

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (vtkIdType id_cell, layer_cells) {
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
      m_LayerAdjacentBoundaryCodes.insert(cell_code->GetValue(id_neigh));
    }
  }
  m_LayerAdjacentBoundaryCodes = m_LayerAdjacentBoundaryCodes - m_BoundaryCodes;

  computeBoundaryLayerVectors();
}

void CreateBoundaryLayer::finalise()
{
  /*
  // set all volume codes
  {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isVolume(id_cell, m_Grid)) {
        cell_code->SetValue(id_cell, m_VolDef.getVC());
      }
    }
  }

  // set m_Grid to modified selected volume + unselected volumes
  {
    MeshPartition volume(m_Grid, true);
    MeshPartition rest(m_RestGrid, true);
    volume.addPartition(rest);
  }
  resetOrientation(m_Grid);
  createIndices(m_Grid);
  */
}

void CreateBoundaryLayer::operate()
{
  prepare();
  writeBoundaryLayerVectors("blayer");

  // ...

  //finalise();
}

