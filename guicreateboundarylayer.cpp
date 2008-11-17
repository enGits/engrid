//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "seedsimpleprismaticlayer.h"
#include "gridsmoother.h"
#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"

void GuiCreateBoundaryLayer::before()
{
  populateBoundaryCodes(ui.listWidget, grid);
};

void GuiCreateBoundaryLayer::operate()
{
  getSelectedItems(ui.listWidget, boundary_codes);
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
          };
        };
      };
      if (bcs.size() >= 2) {
        node_status->SetValue(id_node, 1);
      };
      node_layer->SetValue(id_node, -1);
    };
    foreach (vtkIdType id_cell, layer_cells) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        node_layer->SetValue(pts[i_pts], 0);
      };
    };
  };
  
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
  
  int N_steps = 20;
  seed_layer.setGrid(grid);
  del();
  vol();
  seed_layer.setAllCells();
  seed_layer.setLayerCells(layer_cells);
  seed_layer.setBoundaryCodes(boundary_codes);
  seed_layer();
  seed_layer.getLayerCells(layer_cells);
  
  for (int j = 0; j < N_steps; ++j) {
    cout << "improving prismatic layer -> iteration " << j+1 << "/" << N_steps << endl;
    smooth.setAllCells();
    smooth();
    del.setAllCells();
    del();
    swap();
    vol.setTraceCells(layer_cells);
    vol();
    vol.getTraceCells(layer_cells);
    if (smooth.improvement() < 0.01) break;
  };
  //smooth.setAllCells();
  //smooth();
  
  /*
  del();
  vol();
  cout << "finalising prismatic layer" << endl;
  smooth.setAllCells();
  smooth();
  */
  createIndices(grid);
};

