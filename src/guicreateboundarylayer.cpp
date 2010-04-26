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
#include "removepoints.h"
#include "egvtkobject.h"

GuiCreateBoundaryLayer::GuiCreateBoundaryLayer()
{
  getSet("boundary layer", "number of smoothing iterations", 10, m_NumIterations);

  connect(ui.pushButton_SelectAll_BC, SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(ui.pushButton_ClearAll_BC, SIGNAL(clicked()), this, SLOT(ClearAll_BC()));
}

void GuiCreateBoundaryLayer::before()
{
  ui.checkBoxImprove->setChecked(false);
  l2g_t cells = m_Part.getCells();
  foreach (vtkIdType id_cell, cells) {
    if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
      ui.checkBoxImprove->setChecked(true);
      break;
    }
  }
  populateBoundaryCodes(ui.listWidgetBC);
  populateVolumes(ui.listWidgetVC);
  ui.spinBoxIterations->setValue(m_NumIterations);
  double hr, ha, b;
  getSet("boundary layer", "relative height of boundary layer", 1.0, hr);
  getSet("boundary layer", "absolute height of boundary layer", 1.0, ha);
  getSet("boundary layer", "blending between absolute and relative", 0.0, b);
  {
    int hi = 20*hr;
    hr = 0.05*hi;
    ui.doubleSpinBoxHeight->setValue(hr);
  }
  {
    QString num;
    num.setNum(ha);
    ui.lineEditAbsolute->setText(num);
  }
  {
    int bi = 20*b;
    b = 0.05*bi;
    ui.doubleSpinBoxBlending->setValue(b);
  }
}

void GuiCreateBoundaryLayer::operate()
{
  ///////////////////////////////////////////////////////////////
  // set m_Grid to selected volume
  getSelectedItems(ui.listWidgetBC, m_BoundaryCodes); // fill m_BoundaryCodes with values from listWidgetBC
  QString volume_name = getSelectedVolume(ui.listWidgetVC);
  VolumeDefinition V = GuiMainWindow::pointer()->getVol(volume_name);
  foreach (int bc, m_BoundaryCodes) {
    qDebug()<<"V.getSign("<<bc<<")="<<V.getSign(bc);
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

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach(vtkIdType id_cell, layer_cells) {
      foreach(int i_cell_neighbour, c2c[_cells[id_cell]]) {
          m_LayerAdjacentBoundaryCodes.insert(cell_code->GetValue(cells[i_cell_neighbour]));
      }
  }
  qWarning() << "m_LayerAdjacentBoundaryCodes =" << m_LayerAdjacentBoundaryCodes;
  m_LayerAdjacentBoundaryCodes = m_LayerAdjacentBoundaryCodes - m_BoundaryCodes;

//  m_LayerAdjacentBoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();

  qWarning() << "m_LayerAdjacentBoundaryCodes =" << m_LayerAdjacentBoundaryCodes;
//  EG_BUG;

  cout << "\n\ncreating boundary layer mesh)" << endl;
  
  if (!ui.checkBoxImprove->isChecked()) {
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
  
  if (!ui.checkBoxImprove->isChecked()) {
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
  
  double Hr = ui.doubleSpinBoxHeight->value();
  double Ha = ui.lineEditAbsolute->text().toDouble();
  double bl = ui.doubleSpinBoxBlending->value();
  smooth.setRelativeHeight(Hr);
  smooth.setAbsoluteHeight(Ha);
  smooth.setBlending(bl);
  for (int j = 0; j < ui.spinBoxIterations->value(); ++j) {
    cout << "improving prismatic layer -> iteration " << j+1 << "/" << ui.spinBoxIterations->value() << endl;
    smooth.setAllCells();
    smooth();
    del.setAllCells();
    saveGrid(m_Grid,"del_before");
    del();// does not delete prismatic boundary layer! (->remove points must handle wedges)

    saveGrid(m_Grid,"del_after");
    remove_points();
    cout << "!!!!!!!!!!!!!! REMOVED POINTS: " << remove_points.getNumRemoved() << endl;
//    if(remove_points.getNumRemoved()>0) EG_BUG;

    saveGrid(m_Grid,"after_removepoints");
//    EG_BUG;

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

void GuiCreateBoundaryLayer::SelectAll_BC()
{
  for (int i = 0; i < ui.listWidgetBC->count(); ++i) {
    ui.listWidgetBC->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiCreateBoundaryLayer::ClearAll_BC()
{
  for (int i = 0; i < ui.listWidgetBC->count(); ++i) {
    ui.listWidgetBC->item(i)->setCheckState(Qt::Unchecked);
  }
}
