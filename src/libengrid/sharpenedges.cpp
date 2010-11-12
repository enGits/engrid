//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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

#include "sharpenedges.h"
#include "guimainwindow.h"
#include "laplacesmoother.h"

SharpenEdges::SharpenEdges()
{
  EG_TYPENAME;
  m_PerformGeometricTests             = true;
  m_UseProjectionForSmoothing         = true;
  m_UseNormalCorrectionForSmoothing   = false;
  m_AllowFeatureEdgeSwapping          = false;
  m_RespectFeatureEdgesForDeleteNodes = false;
  m_FeatureAngle                      = deg2rad(180);
  m_EdgeAngle                         = deg2rad(180);
  m_FeatureAngleForDeleteNodes        = deg2rad(20);
  m_NumDelaunaySweeps                 = 5;
  m_NumSmoothSteps                    = 5;
}


void SharpenEdges::assignBCs()
{
  QVector<int> bc_set(m_Grid->GetNumberOfCells(),0);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    bool certain_bc = true;
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      vtkIdType id_node = pts[i];
      for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
        if (cell_code->GetValue(id_cell) != cell_code->GetValue(m_Part.n2cGG(pts[i], j))) {
          certain_bc = false;
          break;
        }
      }
    }
    if (certain_bc) {
      bc_set[id_cell] = cell_code->GetValue(id_cell);
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (bc_set[id_cell] == 0) {
      cell_code->SetValue(id_cell, 9999);
    }
  }
  GuiMainWindow::pointer()->storeSurfaceProjection();
  int N;
  int count = 0;
  double scal_par = 1.0;
  do {
    N = 0;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (bc_set[id_cell] == 0) {
        ++N;
        vec3_t n1 = cellNormal(m_Grid, id_cell);
        int bc = 9999;
        double scal_max = -3;
        vtkIdType id_copy = -1;
        for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
          vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
          vec3_t n2 = cellNormal(m_Grid, id_neigh);
          n1.normalise();
          n2.normalise();
          double scal = n1*n2;
          bool parallel = false;
          if (scal >= scal_par && cell_code->GetValue(id_neigh) < 9999) {
            scal = 1.0;
            parallel = true;
          }
          if (scal > scal_max || parallel) {
            bc = cell_code->GetValue(id_neigh);
            scal_max = scal;
          }
        }
        if (bc < 9999) {
          bc_set[id_cell] = bc;
          cell_code->SetValue(id_cell, bc);
        }
      }
    }
    scal_par *= 0.99;
    cout << scal_par << ',' << N << endl;
    ++count;
  } while (N > 0 && count < 1000);
}


/*
void SharpenEdges::assignBCs()
{
  QVector<int> bc_set(m_Grid->GetNumberOfCells(),0);
  QVector<int> bc_old(m_Grid->GetNumberOfCells(),0);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    bool certain_bc = true;
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    bc_old[id_cell] = cell_code->GetValue(id_cell);
    for (int i = 0; i < N_pts; ++i) {
      for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
        if (cell_code->GetValue(id_cell) != cell_code->GetValue(m_Part.n2cGG(pts[i], j))) {
          certain_bc = false;
          break;
        }
      }
    }
    if (certain_bc) {
      bc_set[id_cell] = cell_code->GetValue(id_cell);
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (bc_set[id_cell] == 0) {
      cell_code->SetValue(id_cell, 9999);
    }
  }
  GuiMainWindow::pointer()->storeSurfaceProjection();

  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (cell_code->GetValue(id_cell) == 9999) {
      cell_code->SetValue(id_cell, bc_old[id_cell]);
    }
  }

  int N;
  int count = 0;
  do {
    N = 0;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      bool ok = false;
      int bc;
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        bc = cell_code->GetValue(m_Part.c2cGG(id_cell, i));
        if (bc == cell_code->GetValue(id_cell)) {
          ok = true;
          break;
        }
      }
      if (!ok) {
        cell_code->SetValue(id_cell, bc);
        ++N;
      }
    }
    ++count;
  } while (N > 0 && count < 1000);
}
*/

void SharpenEdges::operate()
{
  cout << "sharpening edges" << endl;
  assignBCs();
  return;


  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QVector<QSet<int> > n2bc(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    int bc = cell_code->GetValue(id_cell);
    for (int i = 0; i < N_pts; ++i) {
      n2bc[pts[i]].insert(bc);
    }
  }
  int N = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (n2bc[id_node].size() > 1) {
      ++N;
    }
  }
  int total_count = 0;
  int count = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (n2bc[id_node].size() > 1) {
      ++count;
      ++total_count;
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      for (int i_proj_iter = 0; i_proj_iter < 20; ++i_proj_iter) {
        foreach (int bc, n2bc[id_node]) {
          x = GuiMainWindow::pointer()->getSurfProj(bc)->projectFree(x, id_node);
        }
      }
      m_Grid->GetPoints()->SetPoint(id_node, x.data());
      if (count >= 100) {
        cout << total_count << " of " << N << " nodes corrected" << endl;
        count = 0;
      }
    }
  }
  GuiMainWindow::pointer()->storeSurfaceProjection();


/*
  {
    LaplaceSmoother lap;
    lap.setGrid(m_Grid);
    QVector<vtkIdType> cls;
    m_BoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();
    getSurfaceCells(m_BoundaryCodes, cls, m_Grid);
    lap.setCells(cls);
    lap.setNumberOfIterations(5);
    lap.setBoundaryCodes(m_BoundaryCodes);//IMPORTANT: so that unselected nodes become fixed when node types are updated!
    lap.setProjectionOn();
    lap.setNormalCorrectionOff();
    lap.setFreeProjectionForEdgesOn();
    lap();
  }
*/

  UpdateCellIndex(m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
  GuiMainWindow::pointer()->updateActors();
}
