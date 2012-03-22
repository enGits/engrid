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

#include "meshpartition.h"
#include "guimainwindow.h"

#include <vtkKdTreePointLocator.h>
#include <vtkClipDataSet.h>

MeshPartition::MeshPartition()
{
  m_Grid = NULL;
  resetTimeStamps();
}

MeshPartition::MeshPartition(vtkUnstructuredGrid *grid, bool use_all_cells)
{
  resetTimeStamps();
  m_Grid = grid;
  if (use_all_cells) {
    QVector<vtkIdType> cls(grid->GetNumberOfCells());
    for (vtkIdType id_cell = 0; id_cell < cls.size(); ++id_cell) {
      cls[id_cell] = id_cell;
    }
    setCells(cls);
  }
}

MeshPartition::MeshPartition(QString volume_name)
{
  resetTimeStamps();
  setVolume(volume_name);
}

void MeshPartition::resetTimeStamps()
{
  m_CellsStamp = 0;
  m_LCellsStamp = 0;
  m_NodesStamp = 0;
  m_LNodesStamp = 0;
  m_N2NStamp = 0;
  m_N2CStamp = 0;
  m_N2BCStamp = 0;
  m_C2CStamp = 0;
}

void MeshPartition::setVolume(QString volume_name)
{
  m_Grid = GuiMainWindow::pointer()->getGrid();
  resetOrientation(m_Grid);
  VolumeDefinition V = GuiMainWindow::pointer()->getVol(volume_name);
  QList<vtkIdType> cls;
  EG_VTKDCC(vtkIntArray, cell_code,   m_Grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, m_Grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, m_Grid, "cell_voldir");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      int bc = cell_code->GetValue(id_cell);
      cell_voldir->SetValue(id_cell, 0);
      if (V.getSign(bc) != 0) {
        cls.append(id_cell);
        if (V.getSign(bc) == -1) {
          cell_voldir->SetValue(id_cell, 1);
        }
      }
    } else {
      if (cell_code->GetValue(id_cell) == V.getVC()) {
        cls.append(id_cell);
      }
    }
  }
  setCells(cls);
}

void MeshPartition::setVolumeOrientation()
{
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, m_Grid, "cell_voldir");
  foreach (vtkIdType id_cell, m_Cells) {
    if (isSurface(id_cell, m_Grid)) {
      if (cell_curdir->GetValue(id_cell) != cell_voldir->GetValue(id_cell)) {
        reorientateFace(m_Grid, id_cell);
      }
    }
  }
}

void MeshPartition::setOriginalOrientation()
{
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_orgdir, m_Grid, "cell_orgdir");
  foreach (vtkIdType id_cell, m_Cells) {
    if (isSurface(id_cell, m_Grid)) {
      if (cell_curdir->GetValue(id_cell) != cell_orgdir->GetValue(id_cell)) {
        reorientateFace(m_Grid, id_cell);
      }
    }
  }
}

void MeshPartition::setRemainder(const MeshPartition& part)
{
  setGrid(part.getGrid());
  QVector<vtkIdType> rcells;
  getRestCells(m_Grid, part.m_Cells, rcells);
  setCells(rcells);
}

void MeshPartition::extractToVtkGrid(vtkUnstructuredGrid *new_grid)
{
  checkLNodes();
  allocateGrid(new_grid, m_Cells.size(), m_Nodes.size());
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    vec3_t x;
    m_Grid->GetPoints()->GetPoint(m_Nodes[i_nodes], x.data());
    new_grid->GetPoints()->SetPoint(i_nodes, x.data());
    copyNodeData(m_Grid, m_Nodes[i_nodes], new_grid, i_nodes);
  }
  foreach (vtkIdType id_cell, m_Cells) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vtkIdType> new_pts(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      new_pts[i] = m_LNodes[pts[i]];
    }
    vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
    copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
  }
}

void MeshPartition::addPartition(const MeshPartition& part, double tol)
{
  if (m_Grid == part.m_Grid) {
    QVector<bool> cell_marked(m_Grid->GetNumberOfCells(), false);
    foreach (vtkIdType id_cell, m_Cells) {
      cell_marked[id_cell] = true;
    }
    foreach (vtkIdType id_cell, part.m_Cells) {
      cell_marked[id_cell] = true;
    }
    QList<vtkIdType> new_cells;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (cell_marked[id_cell]) {
        new_cells.append(id_cell);
      }
    }
    setCells(new_cells);
  } else {
    if (tol < 0) {
      tol *= -min(getSmallestEdgeLength(), part.getSmallestEdgeLength());
    }
    cout << "  tol=" << tol << endl;
    EG_VTKSP(vtkUnstructuredGrid, new_grid);
    EG_VTKSP(vtkKdTreePointLocator,loc);
    loc->SetDataSet(m_Grid);
    loc->BuildLocator();
    QVector<vtkIdType> pnode2node(part.m_Grid->GetNumberOfPoints());
    vtkIdType N = m_Grid->GetNumberOfPoints();
    Timer T(10);
    for (vtkIdType id_pnode = 0; id_pnode < part.m_Grid->GetNumberOfPoints(); ++id_pnode) {
      vec3_t xp, x;
      part.m_Grid->GetPoint(id_pnode, xp.data());
      vtkIdType id_node = loc->FindClosestPoint(xp.data());
      m_Grid->GetPoint(id_node, x.data());
      if ((x - xp).abs() < tol) {
        pnode2node[id_pnode] = id_node;
      } else {
        pnode2node[id_pnode] = N;
        ++N;
      }
    }
    allocateGrid(new_grid, m_Grid->GetNumberOfCells() + part.m_Cells.size(), N);
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      new_grid->GetPoints()->SetPoint(id_node, x.data());
      copyNodeData(m_Grid, id_node, new_grid, id_node);
    }
    QVector<vtkIdType> part_nodes;
    getNodesFromCells(part.m_Cells, part_nodes, part.m_Grid);
    foreach (vtkIdType id_pnode, part_nodes) {
      vec3_t x;
      part.m_Grid->GetPoint(id_pnode, x.data());
      new_grid->GetPoints()->SetPoint(pnode2node[id_pnode], x.data());
      copyNodeData(part.m_Grid, id_pnode, new_grid, pnode2node[id_pnode]);
    }
    QList<vtkIdType> new_cells;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, pts);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    foreach (vtkIdType id_pcell, part.m_Cells) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = part.m_Grid->GetCellType(id_pcell);
      part.m_Grid->GetCellPoints(id_pcell, N_pts, pts);
      QVector<vtkIdType> new_pts(N_pts);
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = pnode2node[pts[i]];
      }
      vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
      copyCellData(part.m_Grid, id_pcell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    makeCopy(new_grid, m_Grid);
  }
}

double MeshPartition::getSmallestEdgeLength() const
{
  double L = 1e99;
  foreach (vtkIdType id_cell, m_Cells) {
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vec3_t> x(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      m_Grid->GetPoint(pts[i], x[i].data());
    }
    if        (type_cell == VTK_TRIANGLE) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[0]).abs());
    } else if (type_cell == VTK_QUAD) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[3]).abs());
      L = min(L, (x[3] - x[0]).abs());
    } else if (type_cell == VTK_TETRA) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[0]).abs());
      L = min(L, (x[3] - x[0]).abs());
      L = min(L, (x[3] - x[1]).abs());
      L = min(L, (x[3] - x[2]).abs());
    } else if (type_cell == VTK_PYRAMID) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[3]).abs());
      L = min(L, (x[3] - x[0]).abs());
      L = min(L, (x[4] - x[0]).abs());
      L = min(L, (x[4] - x[1]).abs());
      L = min(L, (x[4] - x[2]).abs());
      L = min(L, (x[4] - x[3]).abs());
    } else if (type_cell == VTK_WEDGE) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[0]).abs());
      L = min(L, (x[3] - x[4]).abs());
      L = min(L, (x[4] - x[5]).abs());
      L = min(L, (x[5] - x[3]).abs());
      L = min(L, (x[0] - x[3]).abs());
      L = min(L, (x[1] - x[4]).abs());
      L = min(L, (x[2] - x[5]).abs());
    } else if (type_cell == VTK_HEXAHEDRON) {
      L = min(L, (x[0] - x[1]).abs());
      L = min(L, (x[1] - x[2]).abs());
      L = min(L, (x[2] - x[3]).abs());
      L = min(L, (x[3] - x[0]).abs());
      L = min(L, (x[4] - x[5]).abs());
      L = min(L, (x[5] - x[6]).abs());
      L = min(L, (x[6] - x[7]).abs());
      L = min(L, (x[7] - x[4]).abs());
      L = min(L, (x[0] - x[4]).abs());
      L = min(L, (x[1] - x[5]).abs());
      L = min(L, (x[2] - x[6]).abs());
      L = min(L, (x[3] - x[7]).abs());
    }
  }
  return L;
}

bool MeshPartition::hasNeighNode(vtkIdType id_node, vtkIdType id_neigh)
{
  for (int i = 0; i < n2nGSize(id_node); ++i) {
    if (n2nGG(id_node, i) == id_neigh) {
      return true;
    }
  }
  return false;
}

void MeshPartition::createNodeToBC()
{
  EG_VTKDCC(vtkIntArray, cell_code,   m_Grid, "cell_code");
  m_N2BC.resize(m_Nodes.size());
  foreach (int i_node, m_Nodes) {
    QSet<int> bcs;
    for (int j = 0; j < n2cLSize(i_node); ++j) {
      vtkIdType id_cell = n2cLG(i_node, j);
      if (isSurface(id_cell, m_Grid)) {
        bcs.insert(cell_code->GetValue(n2cLG(i_node, j)));
      }
    }
    m_N2BC[i_node].resize(bcs.size());
    foreach (int bc, bcs) {
      m_N2BC[i_node].append(bc);
    }
  }
}

bool MeshPartition::hasBC(vtkIdType id_node, int bc)
{
  bool found = false;
  for (int j = 0; j < n2bcGSize(id_node); ++j) {
    if (n2bcG(id_node, j) == bc) {
      found == true;
      break;
    }
  }
  return found;
}
