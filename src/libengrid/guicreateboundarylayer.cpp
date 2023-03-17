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
#include "seedsimpleprismaticlayer.h"
#include "gridsmoother.h"
#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"
#include "laplacesmoother.h"
#include "updatedesiredmeshdensity.h"
#include "insertpoints.h"

GuiCreateBoundaryLayer::GuiCreateBoundaryLayer()
{
  getSet("boundary layer", "number of smoothing iterations", 10, m_NumIterations);
  getSet("boundary layer", "remove points", true, m_RemovePoints);

  connect(m_Ui.pushButton_SelectAll_BC, SIGNAL(clicked()), this, SLOT(SelectAll_BC()));
  connect(m_Ui.pushButton_ClearAll_BC, SIGNAL(clicked()), this, SLOT(ClearAll_BC()));
}

void GuiCreateBoundaryLayer::before()
{
  m_Ui.checkBoxImprove->setChecked(false);
  l2g_t cells = m_Part.getCells();
  foreach (vtkIdType id_cell, cells) {
    if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
      m_Ui.checkBoxImprove->setChecked(true);
      break;
    }
  }
  m_Ui.checkBoxRemovePoints->setChecked(m_RemovePoints);
  m_Ui.checkBoxSafeMode->setChecked(false);
  populateBoundaryCodes(m_Ui.listWidgetBC);
  populateVolumes(m_Ui.listWidgetVC);
  m_Ui.spinBoxIterations->setValue(m_NumIterations);
  double hr, ha, b, ds = 1.5, fr = 0.8;
  int num_layers = 0;
  int num_hr, num_nr;
  getSet("boundary layer", "relative height of boundary layer", 0.01, hr);
  getSet("boundary layer", "absolute height of boundary layer", 1.0, ha);
  getSet("boundary layer", "blending between absolute and relative", 0.0, b);
  getSet("boundary layer", "number of layer height relax iterations",  5,  num_hr);
  getSet("boundary layer", "number of normal vector relax iterations", 20, num_nr);
  {
    QString blayer_txt = GuiMainWindow::pointer()->getXmlSection("blayer/global");
    QTextStream s(&blayer_txt);
    if (!s.atEnd()) s >> ha;
    if (!s.atEnd()) s >> hr;
    if (!s.atEnd()) s >> b;
    if (!s.atEnd()) s >> ds;
    if (!s.atEnd()) s >> fr;
    if (!s.atEnd()) s >> num_layers;
    if (!s.atEnd()) s >> num_hr;
    if (!s.atEnd()) s >> num_nr;
  }
  {
    int hi = 2000*hr;
    hr = 0.0005*hi;
    m_Ui.doubleSpinBoxHeight->setValue(hr);
  }
  {
    QString num;
    num.setNum(ha);
    m_Ui.lineEditAbsolute->setText(num);
  }
  {
    if (b > 0.5) {
      m_Ui.checkBoxAbsolute->setChecked(true);
    } else {
      m_Ui.checkBoxAbsolute->setChecked(false);
    }
  }
  {
    m_Ui.doubleSpinBoxStretching->setValue(ds);
    m_Ui.doubleSpinBoxFarRatio->setValue(fr);
  }
  m_Ui.spinBoxHeightIterations->setValue(num_hr);
  m_Ui.spinBoxNormalIterations->setValue(num_nr);
  if (num_layers > 0) {
    QString num;
    num.setNum(num_layers);
    m_Ui.lineEditNumLayers->setText(num);
  } else {
    m_Ui.lineEditNumLayers->setText("not yet available");
  }
}

void GuiCreateBoundaryLayer::reduceSurface()
{
  RemovePoints remove_points;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  remove_points.setMeshPartition(part);
  remove_points.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  remove_points.setUpdatePSPOn();
  remove_points.setThreshold(3);
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), true);

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = false;
      }
    }
  }
  for (int layer = 0; layer < 3; ++layer) {
    QVector<bool> tmp = fix;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!tmp[id_node]) {
        for (int i = 0; i < part.n2nGSize(id_node); ++i) {
          fix[part.n2nGG(id_node, i)] = false;
        }
      }
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = true;
      }
    }
  }

  remove_points.fixNodes(fix);
  remove_points();
}

void GuiCreateBoundaryLayer::updateSurface()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/table").replace("\n", " ");
  int row_count = 0;
  int column_count = 0;
  QVector <VertexMeshDensity>  vmd;
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> row_count >> column_count;
    QVector<int> tmp_bcs;
    GuiMainWindow::pointer()->getAllBoundaryCodes(tmp_bcs);
    if (column_count == tmp_bcs.size() + 3) {
      vmd.fill(VertexMeshDensity(), row_count);
      for (int i = 0; i < row_count; ++i) {
        int row, column;
        QString formula;
        foreach (int bc, tmp_bcs) {
          in >> row >> column >> formula;
          vmd[row].BCmap[bc] = formula.toInt();
        }
        in >> row >> column >> formula;
        vmd[row].type = Str2VertexType(formula);
        in >> row >> column >> formula;
        if (formula == "{{{empty}}}") {
          formula = "";
        }
        vmd[i].setNodes(formula);
        in >> row >> column >> formula;
        vmd[i].density = formula.toDouble();
      }
    } else {
      EG_ERR_RETURN(QObject::tr("Mismatch of number of boundary codes!"));
    }
  }
  UpdateDesiredMeshDensity update;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  update.setMeshPartition(part);
  update.setVertexMeshDensityVector(vmd);
  update.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  update();
}

void GuiCreateBoundaryLayer::insertPoints()
{
  InsertPoints insert_points;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  insert_points.setMeshPartition(part);
  insert_points.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  insert_points();
}

void GuiCreateBoundaryLayer::smoothSurface()
{
  LaplaceSmoother smooth;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  smooth.setMeshPartition(part);
  smooth.setNumberOfIterations(5);
  smooth.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);

  QVector<bool> fix(m_Grid->GetNumberOfPoints(), true);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = false;
      }
    }
  }
  for (int layer = 0; layer < 3; ++layer) {
    QVector<bool> tmp = fix;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!tmp[id_node]) {
        for (int i = 0; i < part.n2nGSize(id_node); ++i) {
          fix[part.n2nGG(id_node, i)] = false;
        }
      }
    }
  }

  smooth.fixNodes(fix);
  smooth();
}

void GuiCreateBoundaryLayer::operate()
{
  if (!GuiMainWindow::pointer()->checkCadInterfaces()) {
    GuiMainWindow::pointer()->storeCadInterfaces();
  }
  ///////////////////////////////////////////////////////////////
  // set m_Grid to selected volume
  getSelectedItems(m_Ui.listWidgetBC, m_BoundaryCodes); // fill m_BoundaryCodes with values from listWidgetBC
  QString volume_name = getSelectedVolume(m_Ui.listWidgetVC);
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

  bool delete_nodes = m_Ui.checkBoxRemovePoints->isChecked();

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach(vtkIdType id_cell, layer_cells) {
    foreach(int i_cell_neighbour, c2c[_cells[id_cell]]) {
      m_LayerAdjacentBoundaryCodes.insert(cell_code->GetValue(cells[i_cell_neighbour]));
    }
  }
  m_LayerAdjacentBoundaryCodes = m_LayerAdjacentBoundaryCodes - m_BoundaryCodes;

  cout << "\n\ncreating boundary layer mesh)" << endl;
  
  if (!m_Ui.checkBoxImprove->isChecked()) {
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
      EG_GET_CELL(id_cell, m_Grid);
      for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
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
  QSet<int> swap_codes = getAllBoundaryCodes(m_Grid);
  swap_codes -= m_LayerAdjacentBoundaryCodes;
  swap.setBoundaryCodes(swap_codes);
  swap.setVerboseOff();

  DeleteTetras del;
  del.setGrid(m_Grid);
  
  if (!m_Ui.checkBoxImprove->isChecked()) {
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
  
  double Hr = m_Ui.doubleSpinBoxHeight->value();
  double Ha = m_Ui.lineEditAbsolute->text().toDouble();
  double bl = 0.0;
  if (m_Ui.checkBoxAbsolute->isChecked()) {
    bl = 1.0;
  }
  smooth.setRelativeHeight(Hr);
  smooth.setAbsoluteHeight(Ha);
  smooth.setBlending(bl);
  smooth.setDesiredStretching(m_Ui.doubleSpinBoxStretching->value());
  smooth.setFarRatio(m_Ui.doubleSpinBoxFarRatio->value());
  smooth.setNumHeightRelaxations(m_Ui.spinBoxHeightIterations->value());
  smooth.setNumNormalRelaxations(m_Ui.spinBoxNormalIterations->value());

  {
    QString blayer_txt = "";
    QTextStream s(&blayer_txt);
    s << Ha << " ";
    s << Hr << " ";
    s << bl << " ";
    s << m_Ui.doubleSpinBoxStretching->value() << " ";
    s << m_Ui.doubleSpinBoxFarRatio->value() << " ";
    s << 0 << " ";
    s << m_Ui.spinBoxHeightIterations->value() << " ";
    s << m_Ui.spinBoxNormalIterations->value() << " ";
    GuiMainWindow::pointer()->setXmlSection("blayer/global", blayer_txt);
  }


  for (int j = 0; j < m_Ui.spinBoxIterations->value(); ++j) {
    cout << "improving prismatic layer -> iteration " << j+1 << "/" << m_Ui.spinBoxIterations->value() << endl;
    if (!m_Ui.checkBoxSafeMode->isChecked()) {
      del.setAllCells();
      del();
    }
    if (delete_nodes || m_Ui.checkBoxSafeMode->isChecked()) {
      smooth.forceNormalCalculation();
    }
    smooth.setAllCells();
    smooth();
    if (delete_nodes) {
      //updateSurface();
      reduceSurface();
      //insertPoints();
    }
    swap();
    smoothSurface();
    swap();
    //vol.setTraceCells(layer_cells);
    if (m_Ui.checkBoxSafeMode->isChecked()) {
      vol();
    }
    //vol.getTraceCells(layer_cells);
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
  for (int i = 0; i < m_Ui.listWidgetBC->count(); ++i) {
    m_Ui.listWidgetBC->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiCreateBoundaryLayer::ClearAll_BC()
{
  for (int i = 0; i < m_Ui.listWidgetBC->count(); ++i) {
    m_Ui.listWidgetBC->item(i)->setCheckState(Qt::Unchecked);
  }
}
