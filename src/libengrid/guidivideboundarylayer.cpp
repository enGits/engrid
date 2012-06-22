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
#include "guidivideboundarylayer.h"
#include "math/linsolve.h"

#include "volumedefinition.h"
#include "guimainwindow.h"

void GuiDivideBoundaryLayer::before()
{
  populateBoundaryCodes(m_Ui.listWidgetBC);
  populateVolumes(m_Ui.listWidgetVC);

  getSet("boundary layer", "first critical angle",  180, m_CritAngle1);
  getSet("boundary layer", "second critical angle", 270, m_CritAngle2);

  QString blayer_txt = GuiMainWindow::pointer()->getXmlSection("blayer/global");
  QTextStream s(&blayer_txt);
  double v;
  if (!s.atEnd()) {
    s >> v;
    QString num;
    num.setNum(v);
    m_Ui.lineEditAbsolute->setText(num);
  }
  if (!s.atEnd()) {
    s >> v; // relative height
    m_Ui.doubleSpinBoxHeight->setValue(v);
  }
  if (!s.atEnd()) {
    s >> v; // blending
    m_Ui.doubleSpinBoxBlending->setValue(v);
  }
  if (!s.atEnd()) {
    s >> v;
    m_Ui.doubleSpinBoxStretching->setValue(v);
  }
  if (!s.atEnd()) {
    s >> v;
    m_Ui.doubleSpinBoxFarRatio->setValue(v);
  }
  if (!s.atEnd()) {
    int v;
    s >> v;
    m_Ui.spinBoxLayers->setValue(v);
  }

  m_Ui.doubleSpinBoxCritAngle1->setValue(m_CritAngle1);
  m_Ui.doubleSpinBoxCritAngle2->setValue(m_CritAngle2);

  m_RestGrid = vtkUnstructuredGrid::New();
}

bool GuiDivideBoundaryLayer::findBoundaryLayer()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");

  l2g_t cells = getPartCells();
  l2l_t c2c   = getPartC2C();

  m_BoundaryCodes.clear();
  m_Pairs.clear();
  m_InsertCell.fill(true,m_Grid->GetNumberOfCells());
  m_NumPrisms = 0;
  m_NumQuads = 0;
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    if (m_Grid->GetCellType(cells[i_cells]) == VTK_WEDGE) {
      ++m_NumPrisms;
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(cells[i_cells],N_pts,pts);
      for (int j = 0; j < 3; ++j) {
        m_Pairs.insert(QPair<vtkIdType,vtkIdType>(pts[j],pts[j+3]));
      }
      for (int j = 2; j < 5; ++j) {
        if (c2c[i_cells][j] != -1) {
          vtkIdType type_ncell = m_Grid->GetCellType(cells[c2c[i_cells][j]]);
          if ((type_ncell != VTK_WEDGE) && (type_ncell != VTK_QUAD)) {
            finalise();
            EG_ERR_RETURN("unable to identify boundary layer");
            return(false);
          }
        } else {
          vec3_t x;
          m_Grid->GetPoint(pts[0],x.data());
          cout << x << endl;
          EG_BUG;
        }
      }
      {
        vtkIdType id_tri = m_Part.c2cLG(i_cells, 0);
        if (m_Grid->GetCellType(id_tri) != VTK_TRIANGLE) {
          EG_BUG;
        }
        m_BoundaryCodes.insert(bc->GetValue(id_tri));
      }
      for (int j = 0; j < 2; ++j) {
        if (c2c[i_cells][j] != -1) {
          vtkIdType type_ncell = m_Grid->GetCellType(cells[c2c[i_cells][j]]);
          if (type_ncell == VTK_WEDGE) {
            finalise();
            EG_ERR_RETURN("the boundary layer seems to have been split already");
            return(false);
          }
        } else {
          EG_BUG;
        }
      }
    }
    if (m_Grid->GetCellType(cells[i_cells]) == VTK_QUAD) {
      ++m_NumQuads;
    }
  }
  if (m_NumPrisms == 0) {
    finalise();
    EG_ERR_RETURN("unable to identify boundary layer");
    return(false);
  }
  m_IsBlayerNode.clear();
  m_IsBlayerNode.fill(false, m_Grid->GetNumberOfPoints());
  
  QPair<vtkIdType,vtkIdType> P;
  foreach (P, m_Pairs) {
    m_IsBlayerNode[P.second] = true;
  }
  
  computeMaxConvexAngles();
  return(true);
}

void GuiDivideBoundaryLayer::computeY1()
{
  /*
  double s1 = 0.01;
  double s2 = 10*m_DesiredStretching;
  while (fabs(s1-s2) > 1e-4) {
    double s = 0.5*(s1+s2);
    for (int i = 2; i < m_Y.size(); ++i) {
      m_Y[i] = m_Y[i-1] + s*(m_Y[i-1]-m_Y[i-2]);
    }
    if (m_Y.last() < 1) {
      s1 = s;
    } else {
      s2 = s;
    }
  }
  */
  double C1 = 0.0;
  double C2 = 2.0;
  while (C2 - C1 > 1e-6) {
    double s = m_DesiredStretching;
    for (int i = 2; i <= m_NumLayers; ++i) {
      m_Y[i] = m_Y[i-1] + s*(m_Y[i-1] - m_Y[i-2]);
      s *= 0.5*(C1 + C2);
    }
    if (m_Y[m_NumLayers] > 1) {
      C2 = 0.5*(C1 + C2);
    } else {
      C1 = 0.5*(C1 + C2);
    }
  }
  m_Y.last() = 1;
}

void GuiDivideBoundaryLayer::computeY2()
{
  double C1 = 0.0;
  double C2 = 100.0;
  double y_target = m_Y[m_NumLayers - 1];
  while (C2 - C1 > 1e-6) {
    double s = 0.5*(C1 + C2);
    for (int i = 2; i < m_NumLayers; ++i) {
      m_Y[i] = m_Y[i-1] + s*(m_Y[i-1] - m_Y[i-2]);
    }
    if (m_Y[m_NumLayers - 1] > y_target) {
      C2 = 0.5*(C1 + C2);
    } else {
      C1 = 0.5*(C1 + C2);
    }
  }
  //m_Y[m_NumLayers - 1] = y_target;
}

void GuiDivideBoundaryLayer::createEdges(vtkUnstructuredGrid *new_grid)
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired" );
  m_Edges.fill(QVector<vtkIdType>(m_NumLayers+1), m_Pairs.size());
  m_Old2Edge.fill(-1, m_Grid->GetNumberOfPoints());
  int N = 0;
  vtkIdType id_new_node = m_Grid->GetNumberOfPoints();
  QPair<vtkIdType,vtkIdType> P;
  double max_step = 0;
  double ymax     = 0;
  double ymin     = 1e99;
  foreach (P, m_Pairs) {
    m_Edges[N][0] = P.first;
    m_Edges[N][m_NumLayers] = P.second;
    m_Old2Edge[P.first] = N;
    m_Old2Edge[P.second] = N;
    
    vec3_t x1,x2;
    m_Grid->GetPoint(P.first, x1.data());
    m_Grid->GetPoint(P.second, x2.data());
    vec3_t n = x2-x1;
    double alpha = GeometryTools::rad2deg(m_MaxConvexAngle[P.first]);
    double h_rel = m_RelativeHeight*cl->GetValue(P.first)/n.abs();
    {
      m_Y.resize(m_NumLayers + 1);
      m_Y[0] = 0;
      m_Y[1] = m_Blending*m_AbsoluteHeight/n.abs() + (1-m_Blending)*h_rel;
      computeY1();
    }
    if (alpha > m_CritAngle1) {
      double blend = min(1.0, (alpha - m_CritAngle1)/(m_CritAngle2 - m_CritAngle1));
      double far_ratio = blend*m_FarRatio + (1-blend)*(1.0 - m_Y[m_NumLayers - 1])*n.abs()/cl->GetValue(P.first);
      m_Y.resize(m_NumLayers + 1);
      m_Y[0] = 0;
      m_Y[1] = m_Blending*m_AbsoluteHeight/n.abs() + (1-m_Blending)*h_rel;
      m_Y[m_NumLayers - 1] = 1.0 - far_ratio*cl->GetValue(P.first)/n.abs();
      m_Y[m_NumLayers] = 1.0;
      computeY2();
    }
    ymin = min(ymin, m_Y[1]*n.abs());
    ymax = max(ymax, m_Y[1]*n.abs());
    for (int i = 1; i < m_NumLayers; ++i) {
      vec3_t x = x1 + m_Y[i]*n;
      max_step = max(max_step,(m_Y[i+1]-m_Y[i])/(m_Y[i]-m_Y[i-1]));
      new_grid->GetPoints()->SetPoint(id_new_node, x.data());
      m_Edges[N][i] = id_new_node;
      ++id_new_node;
    }
    ++N;
  }
  cout << LINE;
  cout << "maximal increment : " << max_step << endl;
  cout << "min(y) : " << ymin << endl;
  cout << "max(y) : " << ymax << endl;
  cout << LINE;
}

void GuiDivideBoundaryLayer::computeMaxConvexAngles()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  m_MaxConvexAngle.fill(0, m_Grid->GetNumberOfPoints());
  EG_FORALL_CELLS (id_cell1, m_Grid) {
    if (isSurface(id_cell1, m_Grid)) {
      if (m_BoundaryCodes.contains(bc->GetValue(id_cell1))) {
        EG_GET_CELL (id_cell1, m_Grid);
        QVector<vec3_t> x(num_pts + 1);
        for (int i = 0; i < num_pts; ++i) {
          m_Grid->GetPoint(pts[i], x[i].data());
        }
        x.last() = x.first();
        vec3_t n1 = cellNormal(m_Grid, id_cell1);
        n1.normalise();
        vec3_t x1 = cellCentre(m_Grid, id_cell1);
        for (int i = 0; i < num_pts; ++i) {
          vec3_t xe = 0.5*(x[i] + x[i+1]);
          vtkIdType id_cell2 = m_Part.c2cGG(id_cell1, i);
          if (m_BoundaryCodes.contains(bc->GetValue(id_cell2))) {
            vec3_t n2 = cellNormal(m_Grid, id_cell2);
            n2.normalise();
            double alpha = angle(n1, n2);
            vec3_t v1 = x1 - xe;
            vec3_t v2 = cellCentre(m_Grid, id_cell2) - xe;
            vec3_t v = v1 + v2;
            v.normalise();
            double sp = v*(n1+n2);
            if (sp > 0) {
              alpha = M_PI + alpha;
            } else {
              alpha = M_PI - alpha;
            }
            vtkIdType p1 = pts[i];
            vtkIdType p2 = pts[0];
            if (i < num_pts - 1) {
              p2 = pts[i+1];
            }
            m_MaxConvexAngle[p1] = max(m_MaxConvexAngle[p1], alpha);
            m_MaxConvexAngle[p2] = max(m_MaxConvexAngle[p2], alpha);
          }
        }
      }
    }
  }
}

void GuiDivideBoundaryLayer::operate()
{
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
  
  m_NumLayers         = m_Ui.spinBoxLayers->value();
  m_RelativeHeight    = m_Ui.doubleSpinBoxHeight->value();
  m_AbsoluteHeight    = m_Ui.lineEditAbsolute->text().toDouble();
  m_Blending          = m_Ui.doubleSpinBoxBlending->value();
  m_DesiredStretching = m_Ui.doubleSpinBoxStretching->value();
  m_FarRatio          = m_Ui.doubleSpinBoxFarRatio->value();
  m_CritAngle1        = m_Ui.doubleSpinBoxCritAngle1->value();
  m_CritAngle2        = m_Ui.doubleSpinBoxCritAngle2->value();
  cout << "dividing boundary layer into " << m_NumLayers << " layers:" << endl;
  if(findBoundaryLayer()) {
    EG_VTKSP(vtkUnstructuredGrid,new_grid);
    allocateGrid(new_grid, m_Grid->GetNumberOfCells() + (m_NumPrisms + m_NumQuads)*(m_NumLayers-1), m_Grid->GetNumberOfPoints() + m_Pairs.size()*(m_NumLayers-1));
    
    
    EG_VTKDCC(vtkIntArray, old_orgdir, m_Grid, "cell_orgdir");
    EG_VTKDCC(vtkIntArray, old_voldir, m_Grid, "cell_voldir");
    EG_VTKDCC(vtkIntArray, old_curdir, m_Grid, "cell_curdir");
    EG_VTKDCC(vtkIntArray, new_orgdir, new_grid, "cell_orgdir");
    EG_VTKDCC(vtkIntArray, new_voldir, new_grid, "cell_voldir");
    EG_VTKDCC(vtkIntArray, new_curdir, new_grid, "cell_curdir");
    
    int orgdir = -99;
    int curdir = -99;
    int voldir = -99;
    
    // copy existing mesh without prisms and adjacent cells
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
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      bool insert_cell = true;
      if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
        insert_cell = false;
      }
      if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
        insert_cell = false;
        if (orgdir != -99 && old_orgdir->GetValue(id_cell) != orgdir) {
          EG_BUG;
        }
        if (voldir != -99 && old_voldir->GetValue(id_cell) != voldir) {
          EG_BUG;
        }
        if (curdir != -99 && old_curdir->GetValue(id_cell) != curdir) {
          EG_BUG;
        }
        orgdir = old_orgdir->GetValue(id_cell);
        voldir = old_voldir->GetValue(id_cell);
        curdir = old_curdir->GetValue(id_cell);
      }
      if (insert_cell) {
        id_new_cell = new_grid->InsertNextCell(m_Grid->GetCellType(id_cell), N_pts, pts);
        copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
      }
    }
    
    // create divided boundary layer
    createEdges(new_grid);
    
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        for (int i = 0; i < m_NumLayers; ++i) {
          vtkIdType p[6];
          p[0] = m_Edges[m_Old2Edge[pts[0]]][i];
          p[1] = m_Edges[m_Old2Edge[pts[1]]][i];
          p[2] = m_Edges[m_Old2Edge[pts[2]]][i];
          p[3] = m_Edges[m_Old2Edge[pts[0]]][i+1];
          p[4] = m_Edges[m_Old2Edge[pts[1]]][i+1];
          p[5] = m_Edges[m_Old2Edge[pts[2]]][i+1];
          id_new_cell = new_grid->InsertNextCell(VTK_WEDGE, 6, p);
          copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
        }
      }
      if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        if ((m_Old2Edge[pts[0]] != -1) && (m_Old2Edge[pts[1]] != -1) && (m_Old2Edge[pts[2]] != -1) && (m_Old2Edge[pts[3]] != -1)) {
          for (int i = 0; i < m_NumLayers; ++i) {
            vtkIdType p[4];
            p[0] = m_Edges[m_Old2Edge[pts[0]]][i];
            p[1] = m_Edges[m_Old2Edge[pts[1]]][i];
            p[2] = m_Edges[m_Old2Edge[pts[1]]][i+1];
            p[3] = m_Edges[m_Old2Edge[pts[0]]][i+1];
            id_new_cell = new_grid->InsertNextCell(VTK_QUAD, 4, p);
            copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
            new_orgdir->SetValue(id_new_cell, orgdir);
            new_voldir->SetValue(id_new_cell, voldir);
            new_curdir->SetValue(id_new_cell, curdir);
          }
        }
      }
    }
    
    makeCopy(new_grid, m_Grid);

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
}

void GuiDivideBoundaryLayer::finalise()
{
  // set m_Grid to modified selected volume + unselected volumes
  {
    MeshPartition volume(m_Grid, true);
    MeshPartition rest(m_RestGrid, true);
    volume.addPartition(rest);
  }
  resetOrientation(m_Grid);
  createIndices(m_Grid);
  m_RestGrid->Delete();
}

