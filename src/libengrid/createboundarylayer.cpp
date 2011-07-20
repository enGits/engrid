// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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

#include "createboundarylayer.h"

#include "guimainwindow.h"
#include "seedsimpleprismaticlayer.h"
#include "gridsmoother.h"
#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"

CreateBoundaryLayer::CreateBoundaryLayer()
{
  getSet("boundary layer", "number of smoothing iterations", 10, m_NumIterations);
  getSet("boundary layer", "remove points", true, m_RemovePoints);
  m_AbsoluteHeight = 0.0;
  m_RelativeHeight = 1.0;
  m_Blending = 0.0;
}

void CreateBoundaryLayer::operate()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  ///////////////////////////////////////////////////////////////
  // set m_Grid to selected volume
  VolumeDefinition V = GuiMainWindow::pointer()->getVol(m_VolumeName);
  foreach (int bc, m_BoundaryCodes) {
    if (V.getSign(bc) == 0) {
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
    rest.extractToVtkGrid(rest_grid);
    makeCopy(vol_grid, m_Grid);
  }
  setAllCells();

  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2l_t  n2c   = getPartN2C();
  l2l_t  c2c   = getPartC2C();
  getSurfaceCells(m_BoundaryCodes, layer_cells, m_Grid);

  bool delete_nodes = m_RemovePoints;

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach(vtkIdType id_cell, layer_cells) {
    foreach(int i_cell_neighbour, c2c[_cells[id_cell]]) {
      m_LayerAdjacentBoundaryCodes.insert(cell_code->GetValue(cells[i_cell_neighbour]));
    }
  }
  m_LayerAdjacentBoundaryCodes = m_LayerAdjacentBoundaryCodes - m_BoundaryCodes;

  cout << "\n\ncreating boundary layer mesh)" << endl;

  {
    EG_VTKDCN(vtkIntArray, node_status, m_Grid, "node_status");
    EG_VTKDCN(vtkIntArray, node_layer,  m_Grid, "node_layer");
    EG_VTKDCC(vtkIntArray, bc,          m_Grid, "cell_code");
    foreach(vtkIdType id_node, nodes) {
      node_status->SetValue(id_node, 0);
      QSet<int> bcs;
      foreach (int i_neigh_cells, n2c[_nodes[id_node]]) {
        vtkIdType id_neigh_cell = cells[i_neigh_cells];
        if (isSurface(id_neigh_cell, m_Grid)) {
          if (m_BoundaryCodes.contains(bc->GetValue(id_neigh_cell))) {
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
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        node_layer->SetValue(pts[i_pts], 0);
      }
    }
  }


  GridSmoother smooth;
  smooth.setGrid(m_Grid);
  smooth.setBoundaryCodes(m_BoundaryCodes);
  smooth.setLayerAdjacentBoundaryCodes(m_LayerAdjacentBoundaryCodes);

  SeedSimplePrismaticLayer seed_layer;

  CreateVolumeMesh vol;
  vol.setGrid(m_Grid);
  SwapTriangles swap;
  swap.setGrid(m_Grid);
  swap.setBoundaryCodes(m_BoundaryCodes);

  RemovePoints remove_points;
  remove_points.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  remove_points.setUpdatePSPOn();

  DeleteTetras del;
  del.setGrid(m_Grid);

  {
    cout << "preparing prismatic layer" << endl;
    seed_layer.setGrid(m_Grid);
    del();
    vol();
    seed_layer.setAllCells();
    seed_layer.setLayerCells(layer_cells);
    seed_layer.setBoundaryCodes(m_BoundaryCodes);
    seed_layer();
    seed_layer.getLayerCells(layer_cells);
  }

  smooth.setRelativeHeight(m_RelativeHeight);
  smooth.setAbsoluteHeight(m_AbsoluteHeight);
  smooth.setBlending(m_Blending);
  for (int j = 0; j < m_NumIterations; ++j) {
    cout << "improving prismatic layer -> iteration " << j+1 << "/" << m_NumIterations << endl;
    smooth.setAllCells();
    smooth();
    del.setAllCells();
    del();// does not delete prismatic boundary layer! (->remove points must handle wedges)

    if(delete_nodes) {
        remove_points();
        qDebug() << "removed points: " << remove_points.getNumRemoved();
    }

    swap();
    vol.setTraceCells(layer_cells);
    vol();
    vol.getTraceCells(layer_cells);
  }

  {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isVolume(id_cell, m_Grid)) {
        cell_code->SetValue(id_cell, V.getVC());
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // set m_Grid to modified selected volume + unselected volumes
  {
    MeshPartition volume(m_Grid, true);
    MeshPartition rest(rest_grid, true);
    volume.addPartition(rest);
  }
  resetOrientation(m_Grid);
  createIndices(m_Grid);
  ///////////////////////////////////////////////////////////////
}

