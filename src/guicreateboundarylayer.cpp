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
#include "guicreateboundarylayer.h"
#include "guimainwindow.h"
#include "seedsimpleprismaticlayer.h"
#include "gridsmoother.h"
#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"

GuiCreateBoundaryLayer::GuiCreateBoundaryLayer()
{
  getSet("boundary layer", "maximal relative error", 0.01, err_max);
  getSet("boundary layer", "maximal number of smoothing iterations", 10, max_iter);
}

void GuiCreateBoundaryLayer::before()
{
  populateBoundaryCodes(ui.listWidgetBC);
  populateVolumes(ui.listWidgetVC);
}

#define DUMP(GRID,NAME) \
{ \
  QString name = QString(NAME) + ".vtu"; \
  setAllCells(); \
  writeCells(GRID, cells, name); \
}

void GuiCreateBoundaryLayer::operate()
{
  getSelectedItems(ui.listWidgetBC, boundary_codes);
  QString volume_name = getSelectedVolume(ui.listWidgetVC);
  VolumeDefinition V = GuiMainWindow::pointer()->getVol(volume_name);
  foreach (int bc, boundary_codes) {
    if (V.getSign(bc) == 0) {
      QString msg;
      msg.setNum(bc);
      msg = "Boundary code " + msg + " is not part of the volume '" + volume_name +"'.";
      EG_ERR_RETURN(msg);
    }
  }

  EG_VTKSP(vtkUnstructuredGrid, rest_grid);
  {
    EG_VTKSP(vtkUnstructuredGrid, vol_grid);
    MeshPartition volume(volume_name);
    MeshPartition rest(grid);
    rest.setRemainder(volume);
    volume.setVolumeOrientation();
    volume.extractToVtkGrid(vol_grid);
    rest.extractToVtkGrid(rest_grid);
    makeCopy(vol_grid, grid);
  }

  DUMP(rest_grid,"reduced");

  setAllCells();
  getSurfaceCells(boundary_codes, layer_cells, grid);

  cout << "\n\ncreating boundary layer mesh)" << endl;
  
  {
    EG_VTKDCN(vtkIntArray, node_status, grid, "node_status");
    EG_VTKDCN(vtkIntArray, node_layer,  grid, "node_layer");
    EG_VTKDCC(vtkIntArray, bc,          grid, "cell_code");
    foreach(vtkIdType id_node, nodes) {
      node_status->SetValue(id_node, 0);
      QSet<int> bcs;
      foreach (int i_neigh_cells, n2c[_nodes[id_node]]) {
        vtkIdType id_neigh_cell = cells[i_neigh_cells];
        if (isSurface(grid, id_neigh_cell)) {
          if (boundary_codes.contains(bc->GetValue(id_neigh_cell))) {
            bcs.insert(bc->GetValue(id_neigh_cell));
          }
        }
      }
      if (bcs.size() >= 2) {
        node_status->SetValue(id_node, 1);
      }
      node_layer->SetValue(id_node, -1);
    }
    foreach (vtkIdType id_cell, layer_cells) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        node_layer->SetValue(pts[i_pts], 0);
      }
    }
  }
  
  cout << "preparing prismatic layer" << endl;
  
  GridSmoother smooth;
  smooth.setGrid(grid);
  smooth.setBoundaryCodes(boundary_codes);
  smooth.prismsOn();
  //smooth.setNumIterations(5);
  
  SeedSimplePrismaticLayer seed_layer; 
  
  CreateVolumeMesh vol;
  vol.setGrid(grid);
  SwapTriangles swap;
  swap.setGrid(grid);
  swap.setBoundaryCodes(boundary_codes);
  DeleteTetras del;
  del.setGrid(grid);
  
  seed_layer.setGrid(grid);
  del();
  vol();
  seed_layer.setAllCells();
  seed_layer.setLayerCells(layer_cells);
  seed_layer.setBoundaryCodes(boundary_codes);
  seed_layer();
  seed_layer.getLayerCells(layer_cells);
  
  for (int j = 0; j < max_iter; ++j) {
    cout << "improving prismatic layer -> iteration " << j+1 << "/" << max_iter << endl;
    smooth.setAllCells();
    smooth();
    del.setAllCells();
    del();
    swap();
    vol.setTraceCells(layer_cells);
    vol();
    vol.getTraceCells(layer_cells);
    if (smooth.improvement() < err_max) break;
    break;
  }
  //smooth.setAllCells();
  //smooth();
  
  /*
  del();
  vol();
  cout << "finalising prismatic layer" << endl;
  smooth.setAllCells();
  smooth();
  */

  {
    MeshPartition volume(grid, true);
    MeshPartition rest(rest_grid, true);
    DUMP(grid,"grid");
    DUMP(rest_grid,"rest_grid");
    volume.addPartition(rest);
  }
  resetOrientation(grid);
  createIndices(grid);
}

