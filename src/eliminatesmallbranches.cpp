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

#include "eliminatesmallbranches.h"
#include "createvolumemesh.h"
#include "guimainwindow.h"

#include <QTime>


EliminateSmallBranches::EliminateSmallBranches()
{
  m_NumLayers = 1;
  m_NumFillLayers = 10;
}

bool EliminateSmallBranches::needsToBeMarked(vtkIdType id_node, int layer)
{
  if (m_IsSurfaceNode[id_node]) {
    return true;
  } else {
    if (layer < m_NumLayers) {
      for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
        if (needsToBeMarked(m_Part.n2nGG(id_node, i), layer+1)) {
          return true;
        }
      }
      return false;
    }
    return false;
  }
}

void EliminateSmallBranches::unmarkNode(vtkIdType id_node, int layer)
{
  if (layer <= m_NumLayers+2) {
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      if (isVolume(id_cell, grid)) {
        m_DeleteCell[id_cell] = false;
      }
    }
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      unmarkNode(m_Part.n2nGG(id_node, i), layer+1);
    }
  }
}

void EliminateSmallBranches::fill(vtkIdType id_cell)
{
  if (!m_DeleteCell[id_cell] && !m_MainVolumeCell[id_cell]) {
    m_MainVolumeCell[id_cell] = true;
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      vtkIdType id_neigh_cell = m_Part.c2cGG(id_cell, i);
      if (id_neigh_cell != -1) {
        fill(id_neigh_cell);
      }
    }
  }
}

void EliminateSmallBranches::fillLayers()
{
  QVector<bool> main_volume_cell = m_MainVolumeCell;
  for (int i_layer = 0; i_layer < m_NumFillLayers; ++i_layer) {
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      if (m_MainVolumeCell[id_cell]) {
        for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
          vtkIdType id_neigh_cell = m_Part.c2cGG(id_cell, i);
          if (id_neigh_cell != -1) {
            if (isVolume(id_neigh_cell, grid)) {
              main_volume_cell[id_neigh_cell] = true;
            }
          }
        }
      }
    }
    m_MainVolumeCell = main_volume_cell;
  }
}

void EliminateSmallBranches::fillCraters()
{
  QVector<bool> is_mainvol_node(grid->GetNumberOfPoints(), false);
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (m_MainVolumeCell[id_cell]) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        is_mainvol_node[pts[i]] = true;
      }
    }
  }
  QVector<bool> is_mainvol_cell(m_MainVolumeCell.size(), true);
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (!m_MainVolumeCell[id_cell]) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        if (!is_mainvol_node[pts[i]]) {
          is_mainvol_cell[id_cell] = false;
          break;
        }
      }
    }
  }
  m_MainVolumeCell = is_mainvol_cell;
}

void EliminateSmallBranches::operate()
{
  QTime time;
  time = QTime::currentTime();
  CreateVolumeMesh vol;
  vol.setGrid(grid);
  vol();
  //vol();
  setAllCells();
  l2g_t  cells = m_Part.getCells();
  g2l_t _cells = m_Part.getLocalCells();
  g2l_t _nodes = m_Part.getLocalNodes();
  l2l_t  n2c   = m_Part.getN2C();
  l2l_t  c2c   = m_Part.getC2C();
  QVector<vtkIdType> faces;
  getAllSurfaceCells(faces, grid);
  m_DeleteCell.fill(false, grid->GetNumberOfCells());
  m_IsSurfaceNode.fill(false, grid->GetNumberOfPoints());
  foreach(vtkIdType id_face, faces) {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_face, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      m_IsSurfaceNode[pts[i]] = true;
    }
    m_DeleteCell[id_face] = true;
  }
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;

  time = QTime::currentTime();
  cout << "marking cells to be removed" << endl;
  foreach(vtkIdType id_cell, cells) {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      if (needsToBeMarked(pts[i])) {
        m_DeleteCell[id_cell] = true;
        break;
      }
    }
  }
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;

  time = QTime::currentTime();
  cout << "'fill' from largest volume" << endl;
  double vol_max = 0;
  vtkIdType id_largest_cell = -1;
  foreach(vtkIdType id_cell, cells) {
    if (isVolume(id_cell, grid)) {
      if (grid->GetCellType(id_cell) != VTK_TETRA) {
        EG_BUG;
      }
      double vol = GeometryTools::cellVA(grid, id_cell, true);
      if (vol > vol_max && !m_DeleteCell[id_cell]) {
        id_largest_cell = id_cell;
        vol_max = vol;
      }
    }
  }
  if (id_largest_cell == -1) {
    EG_BUG;
  }
  m_MainVolumeCell.fill(false, grid->GetNumberOfCells());
  fill(id_largest_cell);
  cout << "ping" << endl;
  fillLayers();
  cout << "pong" << endl;
  fillCraters();
  cout << "pack" << endl;
  for (int i = 0; i < m_MainVolumeCell.size(); ++i) {
    if (m_MainVolumeCell[i]) {
      m_DeleteCell[i] = false;
    } else {
      m_DeleteCell[i] = true;
    }
  }
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;

  time = QTime::currentTime();
  cout << "saving existing boundary faces" << endl;
  foreach (vtkIdType id_face, faces) {
    vtkIdType id_cell = findVolumeCell(grid, id_face, _nodes, cells, _cells, n2c);
    if (id_cell != -1) {
      if (!m_DeleteCell[id_cell]) {
        m_DeleteCell[id_face] = false;
      }
    }
  }
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;

  time = QTime::currentTime();
  cout << "counting new boundary faces" << endl;
  int num_new_faces = 0;
  foreach(vtkIdType id_cell, cells) {
    if (!m_DeleteCell[id_cell]) {
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        if (m_DeleteCell[m_Part.c2cGG(id_cell, i)]) {
          ++num_new_faces;
        }
      }
    }
  }
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;

  time = QTime::currentTime();
  cout << "creating reduced grid" << endl;
  QVector<vtkIdType> old2new(grid->GetNumberOfPoints(), -1);
  vtkIdType num_new_nodes = 0;
  vtkIdType num_new_cells = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (!m_DeleteCell[id_cell]) {
      ++num_new_cells;
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        if (old2new[pts[i]] == -1) {
          old2new[pts[i]] = num_new_nodes;
          ++num_new_nodes;
        }
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, num_new_cells + num_new_faces, num_new_nodes, true);
  EG_VTKDCC(vtkIntArray, bc, new_grid, "cell_code" );
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    if (old2new[id_node] != -1) {
      vec3_t x;
      grid->GetPoint(id_node, x.data());
      new_grid->GetPoints()->SetPoint(old2new[id_node], x.data());
      copyNodeData(grid, id_node, new_grid, old2new[id_node]);
    }
  }
  EG_VTKDCC(vtkIntArray, cell_orgdir, new_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, new_grid, "cell_voldir");
  EG_VTKDCC(vtkIntArray, cell_curdir, new_grid, "cell_curdir");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (!m_DeleteCell[id_cell]) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        pts[i] = old2new[pts[i]];
        if (pts[i] == -1) {
          EG_BUG;
        }
      }
      vtkIdType id_new_cell = new_grid->InsertNextCell(grid->GetCellType(id_cell), N_pts, pts);
      copyCellData(grid, id_cell, new_grid, id_new_cell);
      if (grid->GetCellType(id_cell) == VTK_TETRA) {
        for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
          if (m_DeleteCell[m_Part.c2cGG(id_cell, i)]) {
            vtkIdType face_pts[3];
            if (i == 0) {
              face_pts[0] = pts[2];
              face_pts[1] = pts[1];
              face_pts[2] = pts[0];
            } else if (i == 1) {
              face_pts[0] = pts[1];
              face_pts[1] = pts[3];
              face_pts[2] = pts[0];
            } else if (i == 2) {
              face_pts[0] = pts[3];
              face_pts[1] = pts[2];
              face_pts[2] = pts[0];
            } else if (i == 3) {
              face_pts[0] = pts[2];
              face_pts[1] = pts[3];
              face_pts[2] = pts[1];
            }
            vtkIdType id_new_cell = new_grid->InsertNextCell(VTK_TRIANGLE, 3, face_pts);
            copyCellData(grid, id_cell, new_grid, id_new_cell);
            bc->SetValue(id_new_cell, 99);
            cell_orgdir->SetValue(id_new_cell, 0);
            cell_voldir->SetValue(id_new_cell, 0);
            cell_curdir->SetValue(id_new_cell, 0);
          }
        }
      }
    }
  }
  makeCopy(new_grid, grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
  cout << time.secsTo(QTime::currentTime()) << " seconds" << endl;
}
