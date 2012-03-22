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
#include "polymesh.h"


// ==========
//   face_t
// ==========

PolyMesh::face_t::face_t(int N, int o, int n, int b)
{
  node.resize(N);
  owner = o;
  neighbour = n;
  bc = b;
}

int PolyMesh::face_t::operator[](int i)
{
  while (i < 0) {
    i += node.size();
  }
  while (i >= node.size()) {
    i -= node.size();
  }
  return node[i];
}

bool PolyMesh::face_t::operator<(const face_t &F) const
{
  bool less = false;
  if (bc < F.bc) {
    less = true;
  } else if (bc == F.bc) {
    if (owner < F.owner) {
      less = true;
    } else if (owner == F.owner) {
      if (neighbour < F.neighbour) {
        less = true;
      }
    }
  }
  return less;
}

// ==========
//   node_t
// ==========

PolyMesh::node_t::node_t(vtkIdType id1)
{
  id.resize(1);
  id[0] = id1;
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2)
{
  id.resize(2);
  id[0] = id1;
  id[1] = id2;
  qSort(id);
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3)
{
  id.resize(3);
  id[0] = id1;
  id[1] = id2;
  id[2] = id3;
  qSort(id);
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType id4)
{
  id.resize(4);
  id[0] = id1;
  id[1] = id2;
  id[2] = id3;
  id[3] = id4;
  qSort(id);
}

bool PolyMesh::node_t::operator<(const PolyMesh::node_t &N) const
{
  int num = min(id.size(), N.id.size());
  for (int i = 0; i < num; ++i) {
    if (id[i] < N.id[i]) {
      return true;
    }
    if (id[i] > N.id[i]) {
      return false;
    }
  }
  if (id.size() < N.id.size()) {
    return true;
  }
  return false;
}

bool PolyMesh::node_t::operator>(const PolyMesh::node_t &N) const
{
  int num = min(id.size(), N.id.size());
  for (int i = 0; i < num; ++i) {
    if (id[i] > N.id[i]) {
      return true;
    }
    if (id[i] < N.id[i]) {
      return false;
    }
  }
  if (id.size() > N.id.size()) {
    return true;
  }
  return false;
}

bool PolyMesh::node_t::operator==(const PolyMesh::node_t &N) const
{
  if (id.size() != N.id.size()) {
    return false;
  }
  for (int i = 0; i < id.size(); ++i) {
    if (id[i] != N.id[i]) {
      return false;
    }
  }
  return true;
}





// ============
//   PolyMesh
// ============

PolyMesh::PolyMesh(vtkUnstructuredGrid *grid, bool dual_mesh)
{
  if (!dual_mesh) {
    EG_BUG;
  }
  m_Grid = grid;
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
  findDualCells();
}

void PolyMesh::getFacesOfEdgeInsideTetra(vtkIdType id_cell, vtkIdType id_node1, vtkIdType id_node2, int &face1, int &face2)
{
  if (m_Grid->GetCellType(id_cell) != VTK_TETRA) {
    EG_BUG;
  }
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell, num_pts, pts);
  face1 = -1;
  face2 = -1;

  if        (id_node1 == pts[0]) {
    if      (id_node2 == pts[1]) { face1 = 0; face2 = 1; }
    else if (id_node2 == pts[2]) { face1 = 2; face2 = 0; }
    else if (id_node2 == pts[3]) { face1 = 1; face2 = 2; }

  } else if (id_node1 == pts[1]) {
    if      (id_node2 == pts[0]) { face1 = 1; face2 = 0; }
    else if (id_node2 == pts[2]) { face1 = 0; face2 = 3; }
    else if (id_node2 == pts[3]) { face1 = 3; face2 = 1; }

  } else if (id_node1 == pts[2]) {
    if      (id_node2 == pts[0]) { face1 = 0; face2 = 2; }
    else if (id_node2 == pts[1]) { face1 = 3; face2 = 0; }
    else if (id_node2 == pts[3]) { face1 = 2; face2 = 3; }

  } else if (id_node1 == pts[3]) {
    if      (id_node2 == pts[0]) { face1 = 2; face2 = 1; }
    else if (id_node2 == pts[1]) { face1 = 1; face2 = 3; }
    else if (id_node2 == pts[2]) { face1 = 3; face2 = 2; }
  }

  if (face1 == -1 || face2 == -1) {
    EG_BUG;
  }
}

void PolyMesh::getSortedEdgeCells(vtkIdType id_node1, vtkIdType id_node2, QList<vtkIdType> &cells)
{
  /*
  QSet<vtkIdType> tetras_node1;
  QSet<vtkIdType> tetras_node2;
  for (int j = 0; j < m_Part.n2bcGSize(id_node1); ++j) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, j);
    if (isVolume(id_cell, m_Grid)) {
      if (m_Grid->GetCellType(id_cell) != VTK_TETRA) {
        EG_BUG;
      }
      tetras_node1.insert(id_cell);
    }
  }
  for (int j = 0; j < m_Part.n2bcGSize(id_node2); ++j) {
    vtkIdType id_cell = m_Part.n2cGG(id_node2, j);
    if (isVolume(id_cell, m_Grid)) {
      if (m_Grid->GetCellType(id_cell) != VTK_TETRA) {
        EG_BUG;
      }
      tetras_node2.insert(id_cell);
    }
  }
  QSet<vtkIdType> tetras_raw = tetras_node1.intersect(tetras_node2);
  QList<vtkIdType> tetras_sorted;
  */

  cells.clear();
  vtkIdType id_start = -1;
  for (int i = 0; i < m_Part.n2bcGSize(id_node1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, i);
    if (isVolume(id_cell, m_Grid)) {
      if (m_Grid->GetCellType(id_cell) != VTK_TETRA) {
        EG_BUG;
      }
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_cell, num_pts, pts);
      for (int j = 0; j < num_pts; ++j) {
        if (pts[i] == id_node2) {
          id_start = id_cell;
          break;
        }
      }
    }
  }
  if (id_start == -1) {
    EG_BUG;
  }
  vtkIdType id_cell = -1;
  do {
  } while (id_cell != id_start && !isSurface(id_cell, m_Grid));

  //Hier gehts weiter !!!
}

void PolyMesh::findDualCells()
{
  m_Cell2PCell.fill(m_Grid->GetNumberOfCells(), -1);
  m_Node2PCell.fill(m_Grid->GetNumberOfPoints(), -1);
  // ...
}

void PolyMesh::countNodesAndFaces()
{
  int num_poly_cells = 0;
  int num_faces = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (!isHexCoreNode(id_node)) {
      int num_vol = 0;
      int num_surf = 0;
      bool tetras_only = true;
      for (int j = 0; j < m_Part.n2cGSize(id_node); ++j) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, j);
        if (isVolume(id_cell, m_Grid)) {
          ++num_vol;
          if (m_Grid->GetCellType(id_cell) != VTK_TETRA) {
            tetras_only = false;
          }
        } else {
          ++num_surf;
        }
      }
      if (tetras_only) {
        ++num_poly_cells;
        num_faces += num_vol + num_surf;
      }
    }
  }
}

