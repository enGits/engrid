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

#include "meshpartition.h"
#include "guimainwindow.h"

#include <vtkKdTreePointLocator.h>
#include <vtkClipDataSet.h>
#include <vtkTriangleFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkSTLWriter.h>

#include <vtkEgNormalExtrusion.h>

MeshPartition::MeshPartition()
{
  m_Grid = NULL;
  m_TrackGrid = false;
  resetTimeStamps();
}

MeshPartition::MeshPartition(vtkUnstructuredGrid *grid, bool use_all_cells)
{
  setGrid(grid, use_all_cells);
}

MeshPartition::MeshPartition(QString volume_name)
{
  m_TrackGrid = false;
  resetTimeStamps();
  setVolume(volume_name);
}

void MeshPartition::setGrid(vtkUnstructuredGrid *grid, bool use_all_cells)
{
  m_TrackGrid = false;
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

void MeshPartition::trackGrid(vtkUnstructuredGrid *grid)
{
  setGrid(grid);
  setAllCells();
  m_GridMTime = m_Grid->GetMTime();
  m_TrackGrid = true;
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
    /*
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vtkIdType> new_pts(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      new_pts[i] = m_LNodes[pts[i]];
    }
    // update for polyhedral cells here
    vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
    */
    vtkIdType id_new_cell = copyCell(m_Grid, id_cell, new_grid, m_LNodes);
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
    tol = max(tol, 1e-30);
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
      vtkIdType id_new_cell = copyCell(m_Grid, id_cell, new_grid);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    foreach (vtkIdType id_pcell, part.m_Cells) {
      vtkIdType id_new_cell = copyCell(part.m_Grid, id_pcell, new_grid, pnode2node);
      copyCellData(part.m_Grid, id_pcell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    makeCopy(new_grid, m_Grid);
  }
}

void MeshPartition::concatenatePartition(const MeshPartition& part)
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
    EG_VTKSP(vtkUnstructuredGrid, new_grid);
    QVector<vtkIdType> pnode2node(part.m_Grid->GetNumberOfPoints());
    vtkIdType N = m_Grid->GetNumberOfPoints();
    for (vtkIdType id_pnode = 0; id_pnode < part.m_Grid->GetNumberOfPoints(); ++id_pnode) {
      pnode2node[id_pnode] = N;
      ++N;
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
      vtkIdType id_new_cell = copyCell(m_Grid, id_cell, new_grid);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
      new_cells.append(id_new_cell);
    }
    foreach (vtkIdType id_pcell, part.m_Cells) {
      vtkIdType id_new_cell = copyCell(part.m_Grid, id_pcell, new_grid, pnode2node);
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
  for (int i_node = 0; i_node < m_Nodes.size(); ++i_node) {
    QSet<int> bcs;
    for (int j = 0; j < n2cLSize(i_node); ++j) {
      vtkIdType id_cell = n2cLG(i_node, j);
      if (isSurface(id_cell, m_Grid)) {
        int bc = cell_code->GetValue(n2cLG(i_node, j));
        if (bc != 0) {
          bcs.insert(cell_code->GetValue(n2cLG(i_node, j)));
        }
      }
    }
    m_N2BC[i_node].resize(bcs.size());
    int i = 0;
    foreach (int bc, bcs) {
      m_N2BC[i_node][i++] = bc;
    }
  }
}

bool MeshPartition::hasBC(vtkIdType id_node, int bc)
{
  bool found = false;
  for (int j = 0; j < n2bcGSize(id_node); ++j) {
    if (n2bcG(id_node, j) == bc) {
      found = true;
      break;
    }
  }
  return found;
}

vtkIdType MeshPartition::getVolumeCell(vtkIdType id_face)
{
  checkLNodes();
  checkLCells();
  checkN2C();
  return findVolumeCell(m_Grid, id_face, m_LNodes, m_Cells, m_LCells, m_N2C);
}

vec3_t MeshPartition::globalNormal(vtkIdType id_node)
{
  vec3_t normal(0,0,0);
  for (int i = 0; i < n2cGSize(id_node); ++i) {
    vtkIdType id_cell = n2cGG(id_node, i);
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t a, b, c;
      for (int j = 0; j < N_pts; ++j) {
        if (pts[j] == id_node) {
          m_Grid->GetPoint(pts[j], a.data());
          if (j > 0) {
            m_Grid->GetPoint(pts[j-1], b.data());
          } else {
            m_Grid->GetPoint(pts[N_pts-1], b.data());
          }
          if (j < N_pts - 1) {
            m_Grid->GetPoint(pts[j+1], c.data());
          } else {
            m_Grid->GetPoint(pts[0], c.data());
          }
        }
      }
      vec3_t u = b - a;
      vec3_t v = c - a;
      double alpha = GeometryTools::angle(u, v);
      vec3_t n = u.cross(v);
      n.normalise();
      if (checkVector(n)) {
        normal -= alpha*n;
      }
    }
  }
  normal.normalise();
  return normal;
}

double MeshPartition::getAverageSurfaceEdgeLength(vtkIdType id_node)
{
  QSet<vtkIdType> surface_neighbours;
  for (int i = 0; i < n2cGSize(id_node); ++i) {
    vtkIdType id_cell = n2cGG(id_node, i);
    if (isSurface(id_cell, m_Grid)) {
      EG_GET_CELL_AND_TYPE(id_cell, m_Grid);
      for (int j = 0; j < num_pts; ++j) {
        if (pts[j] != id_node) {
          surface_neighbours.insert(pts[j]);
        }
      }
    }
  }
  double L = 0;
  if (surface_neighbours.size() > 0) {
    vec3_t x, xn;
    m_Grid->GetPoint(id_node, x.data());
    foreach (vtkIdType id_neigh, surface_neighbours) {
      m_Grid->GetPoint(id_neigh, xn.data());
      L += (x - xn).abs();
    }
    L /= surface_neighbours.size();
  }
  return L;
}

void MeshPartition::computeMinAndMaxSurfaceStencilEdgeLengths(vtkIdType id_node, double &l_min, double &l_max)
{
  l_min = 1e99;
  l_max = 0;
  for (int i = 0; i < n2cGSize(id_node); ++i) {
    vtkIdType id_cell = n2cGG(id_node, i);
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType *pts, num_pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      vec3_t x1, x2;
      m_Grid->GetPoint(pts[num_pts - 1], x1.data());
      for (int j = 0; j < num_pts; ++j) {
        m_Grid->GetPoint(pts[j], x2.data());
        double L = (x1 - x2).abs();
        l_min = min(L, l_min);
        l_max = max(L, l_max);
        x1 = x2;
      }
    }
  }
}

double MeshPartition::getMinSurfaceStencilEdgeLength(vtkIdType id_node)
{
  double l_min, l_max;
  computeMinAndMaxSurfaceStencilEdgeLengths(id_node, l_min, l_max);
  return l_min;
}

double MeshPartition::getMaxSurfaceStencilEdgeLength(vtkIdType id_node)
{
  double l_min, l_max;
  computeMinAndMaxSurfaceStencilEdgeLengths(id_node, l_min, l_max);
  return l_max;
}

int MeshPartition::getNumberOfFeatureNeighbours(vtkIdType id_node)
{
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  int N = 0;
  for (int i = 0; i < n2nGSize(id_node); ++i) {
    char type = node_type->GetValue(n2nGG(id_node, i));
    if (type == EG_FEATURE_EDGE_VERTEX || type == EG_FEATURE_CORNER_VERTEX) {
      ++N;
    }
  }
  return N;
}

int MeshPartition::getEdgeType(vtkIdType id_node1, vtkIdType id_node2)
{
  QList <vtkIdType> edge_faces;
  getEdgeFaces(id_node1, id_node2, edge_faces);
  int edge = 0;
  if (edge_faces.size() == 1) {
    edge = 2;
  } else if (edge_faces.size() >= 2) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    QSet<int> bcs;
    foreach (vtkIdType id_face, edge_faces) {
      bcs.insert(cell_code->GetValue(id_face));
    }
    if (bcs.size() > 1) {
      edge = 2;
    } else {
      edge = 1;
    }
  }
  return edge;
}

int MeshPartition::computeTopoDistance(vtkIdType id_node1, vtkIdType id_node2, int max_dist, int restriction_type)
{
  int dist = 0;
  QSet<vtkIdType> candidates;
  candidates.insert(id_node1);
  while (dist < max_dist && !candidates.contains(id_node2)) {
    foreach (vtkIdType id_node, candidates) {
      for (int i = 0; i < n2nGSize(id_node); ++i) {
        vtkIdType id_neigh = n2nGG(id_node,i);
        if (!candidates.contains(id_neigh)) {
          bool insert = true;
          if (restriction_type > 0) {
            insert = getEdgeType(id_node, id_neigh) >= restriction_type;
          }
          if (insert) {
            candidates.insert(id_neigh);
          }
        }
      }
    }
    ++dist;
  }
  return dist;
}

void MeshPartition::getCommonNodes(vtkIdType id_cell1, vtkIdType id_cell2, QVector<vtkIdType> &common_nodes)
{
  common_nodes.clear();
  QSet<vtkIdType> nodes1, nodes2;
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell1, num_pts, pts);
  for (int i = 0; i < num_pts; ++i) {
    nodes1.insert(pts[i]);
  }
  m_Grid->GetCellPoints(id_cell2, num_pts, pts);
  for (int i = 0; i < num_pts; ++i) {
    nodes2.insert(pts[i]);
  }
  nodes1.intersect(nodes2);
  common_nodes.resize(nodes1.size());
  qCopy(nodes1.begin(), nodes1.end(), common_nodes.begin());
}

bool MeshPartition::isFeatureEdge(vtkIdType id_node1, vtkIdType id_node2, double feature_angle)
{
  bool is_feature_edge = false;
  QVector<int> bcs;
  QVector<vtkIdType> nodes(2);
  nodes[0] = id_node1;
  nodes[1] = id_node2;
  getCommonBcs(nodes, bcs);
  if (bcs.size() == 1) {
    QList<vtkIdType> faces;
    getEdgeFaces(id_node1, id_node2, faces);
    if (faces.size() == 2) {
      vec3_t n1 = cellNormal(m_Grid, faces[0]);
      vec3_t n2 = cellNormal(m_Grid, faces[1]);
      if (GeometryTools::angle(n1, n2) > feature_angle) {
        EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
        CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(bcs.first());
        vec3_t x1, x2;
        m_Grid->GetPoint(id_node1, x1.data());
        m_Grid->GetPoint(id_node2, x2.data());
        vec3_t x = 0.5*(x1 + x2);
        double h = min(cl->GetValue(id_node1), cl->GetValue(id_node2));
        vec3_t xs = cad->snapToEdge(x);
        if ((xs - x).abs() < 0.2*h) {
          is_feature_edge = true;
        }
      }
    }
  }
  return is_feature_edge;
}

bool MeshPartition::isConvexNode(vtkIdType id_node)
{
  int n_faces = n2cGSize(id_node);
  if (n_faces <= 1) {
    return false;
  }

  vec3_t x1, cell_centers(0,0,0), cell_normals(0,0,0);
  m_Grid->GetPoint(id_node, x1.data());
  for (int i = 0; i < n_faces; ++i) {
    vtkIdType id_cell = n2cGG(id_node, i);
    cell_centers += cellCentre(m_Grid, id_cell);
    cell_normals += cellNormal(m_Grid, id_cell);
  }
  cell_centers *= 1.0/n_faces;
  if ((x1 - cell_centers)*cell_normals.normalise() > 0) {
    return true;
  }
  return false;
}


bool MeshPartition::isConvexNode(vtkIdType id_node, QVector<int> bl_codes)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  int n_faces = n2cGSize(id_node);
  if (n_faces <= 1) {
    return false;
  }

  vec3_t x1, cell_centers(0,0,0), cell_normals(0,0,0);
  m_Grid->GetPoint(id_node, x1.data());
  int n_used = 0;
  for (int i = 0; i < n_faces; ++i) {
    vtkIdType id_cell = n2cGG(id_node, i);
    if (bl_codes.contains(cell_code->GetValue(id_cell))) {
      cell_centers += cellCentre(m_Grid, id_cell);
      cell_normals += cellNormal(m_Grid, id_cell);
      ++n_used;
    }
  }
  if (n_used > 0) {
    cell_centers *= 1.0/n_used;
    if ((x1 - cell_centers)*cell_normals.normalise() > 0) {
      return true;
    }
  }
  return false;
}

void MeshPartition::calcPlanarSurfaceMetrics(double &Dh, double &A, double &P, vec3_t &x, vec3_t &n)
{
  A = 0;
  n = vec3_t(0, 0, 0);
  x = vec3_t(0, 0, 0);
  P = 0;
  foreach (vtkIdType id_cell, m_Cells) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      QVector<vec3_t> x_pt(num_pts + 1);
      for (int i = 0; i < num_pts; ++i) {
        m_Grid->GetPoint(pts[i], x_pt[i].data());
      }
      x_pt[num_pts] = x_pt[0];
      double a = cellVA(m_Grid, id_cell);
      A += a;
      x += a*cellCentre(m_Grid, id_cell);
      n += cellNormal(m_Grid, id_cell);
      for (int i = 0; i < num_pts; ++i) {
        if (c2cGG(id_cell, i) == -1) {
          P += (x_pt[i+1] - x_pt[i]).abs();
        }
      }
    }
  }
  n.normalise();
  x *= 1.0/A;
  Dh = 4*A/P;
}

bool MeshPartition::isPlanar(double tolerance_angle)
{
  foreach (vtkIdType id_cell, m_Cells) {
    if (!isSurface(id_cell, m_Grid)) {
      return false;
    }
  }
  double Dh, A, P;
  vec3_t x, n0;
  calcPlanarSurfaceMetrics(Dh, A, P, n0, x);
  foreach (vtkIdType id_cell, m_Cells) {
    vec3_t n = cellNormal(m_Grid, id_cell);
    if (angle(n, n0) > tolerance_angle) {
      return false;
    }
  }
  return true;
}

void MeshPartition::setBC(int bc)
{
  QList<vtkIdType> cls;
  EG_VTKDCC(vtkIntArray, cell_code,   m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (cell_code->GetValue(id_cell) == bc) {
      cls.append(id_cell);
    }
  }
  setCells(cls);
}

void MeshPartition::duplicate()
{
  vtkIdType old_num_cells = m_Grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  extractToVtkGrid(new_grid);
  MeshPartition new_part(new_grid, true);
  concatenatePartition(new_part);
  QList<vtkIdType> cells;
  for (vtkIdType id_cell = old_num_cells; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    cells << id_cell;
  }
  setCells(cells);
}

void MeshPartition::scale(double factor, vec3_t centre)
{
  checkNodes();
  foreach (vtkIdType id_node, m_Nodes) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    x -= centre;
    x *= factor;
    x += centre;
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
}

void MeshPartition::translate(vec3_t v)
{
  checkNodes();
  foreach (vtkIdType id_node, m_Nodes) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    x += v;
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
}

void MeshPartition::extrude(vec3_t dir, QList<double> h, BoundaryCondition extrude_bc,
                            bool force_bottom_bc, bool force_side_bc, bool force_top_bc,
                            BoundaryCondition bottom_bc, BoundaryCondition side_bc, BoundaryCondition top_bc)
{
  checkNodes();
  if (onlySurfaceCells()) {
    QList<vtkIdType> new_cells;
    foreach (vtkIdType id_cell, m_Cells) {
      new_cells << id_cell;
    }
    vtkIdType old_num_cells = m_Grid->GetNumberOfCells();
    extrude_bc = GuiMainWindow::pointer()->getBC(extrude_bc);
    EG_VTKDCC(vtkIntArray, cell_code,   m_Grid, "cell_code");
    foreach (vtkIdType id_cell, m_Cells) {
      cell_code->SetValue(id_cell, extrude_bc.getCode());
    }
    QSet<int> bcs;
    bcs.insert(extrude_bc.getCode());
    EG_VTKSP(vtkEgNormalExtrusion, extr);
    QVector<double> y(h.size() + 1);
    {
      y[0] = 0;
      int i = 1;
      foreach (double dy, h) {
        y[i] = y[i-1] + dy;
        ++i;
      }
    }
    extr->SetLayers(y);
    extr->SetBoundaryCodes(bcs);
    if (force_bottom_bc) extr->SetCustomBottomBc (bottom_bc.getCode());
    if (force_side_bc)   extr->SetCustomSideBc   (side_bc.getCode());
    if (force_top_bc)    extr->SetCustomTopBc    (top_bc.getCode());
    dir.normalise();
    extr->SetNormal(dir);
    extr->SetFixed();
    EG_VTKSP(vtkUnstructuredGrid,ug);
    makeCopy(m_Grid, ug);
    extr->SetInputData(ug);
    extr->Update();
    makeCopy(extr->GetOutput(), m_Grid);
    setBC(extrude_bc.getCode());
  }
}

bool MeshPartition::onlySurfaceCells()
{
  foreach (vtkIdType id_cell, m_Cells) {
    if (!isSurface(id_cell, m_Grid)) {
      return false;
    }
  }
  return true;
}

void MeshPartition::resetBC(QString bc_name, QString bc_type)
{
  BoundaryCondition bc = GuiMainWindow::pointer()->getBC(BoundaryCondition(bc_name, bc_type));
  EG_VTKDCC(vtkIntArray, cell_code,   m_Grid, "cell_code");
  foreach (vtkIdType id_cell, m_Cells) {
    cell_code->SetValue(id_cell, bc.getCode());
  }
}

void MeshPartition::writeSTL(QString file_name)
{
  if (file_name.right(4).toLower() != ".stl") {
    file_name += ".stl";
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  extractToVtkGrid(grid);
  EG_VTKSP(vtkGeometryFilter, geometry);
  geometry->SetInputData(grid);
  EG_VTKSP(vtkTriangleFilter, triangle);
  triangle->SetInputConnection(geometry->GetOutputPort());
  EG_VTKSP(vtkSTLWriter, stl);
  stl->SetInputConnection(  triangle->GetOutputPort());
  stl->SetFileName(qPrintable(file_name));
  stl->SetFileTypeToBinary();
  stl->Write();
}
