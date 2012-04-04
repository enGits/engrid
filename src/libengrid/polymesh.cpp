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

PolyMesh::face_t::face_t(int N, int o, int n, vec3_t rv, int b)
{
  node.resize(N);
  owner = o;
  neighbour = n;
  bc = b;
  ref_vec = rv;
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

bool PolyMesh::face_t::operator==(const face_t &F) const
{
  bool equal = false;
  if (bc == F.bc) {
    if (owner == F.owner) {
      if (neighbour == F.neighbour) {
        equal = true;
      }
    }
  }
  return equal;
}


// ==========
//   node_t
// ==========

PolyMesh::node_t::node_t(const QVector<vtkIdType> &ids)
{
  id = ids;
  qSort(id);
}

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
  findPolyCells();
  createNodesAndFaces();
  checkFaceOrientation();

  /*
  for (int i = 0; i < numCells(); ++i) {
    vec3_t n_total(0,0,0);
    for (int j = 0; j < m_Faces.size(); ++j) {
      vec3_t n = faceNormal(j);
      bool found = false;
      if (m_Faces[j].owner == i) {
        n_total += n;
        found = true;
      }
      if (m_Faces[j].neighbour == i) {
        n_total -= n;
        found = true;
      }
      if (found) {
        cout << "break" << endl;
      }
    }
    cout << i << ", " << n_total << endl;
  }
  */
}

void PolyMesh::getFacesOfEdgeInsideCell(vtkIdType id_cell, vtkIdType id_node1, vtkIdType id_node2, int &face1, int &face2)
{
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell, num_pts, pts);
  face1 = -1;
  face2 = -1;

  if (m_Grid->GetCellType(id_cell) == VTK_TETRA) {
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
  } else if (m_Grid->GetCellType(id_cell) == VTK_PYRAMID) {
    if        (id_node1 == pts[0]) {
      if      (id_node2 == pts[1]) { face1 = 0; face2 = 1; }
      else if (id_node2 == pts[3]) { face1 = 4; face2 = 0; }
      else if (id_node2 == pts[4]) { face1 = 1; face2 = 4; }
    } else if (id_node1 == pts[1]) {
      if      (id_node2 == pts[0]) { face1 = 1; face2 = 0; }
      else if (id_node2 == pts[2]) { face1 = 0; face2 = 2; }
      else if (id_node2 == pts[4]) { face1 = 2; face2 = 1; }
    } else if (id_node1 == pts[2]) {
      if      (id_node2 == pts[1]) { face1 = 2; face2 = 0; }
      else if (id_node2 == pts[3]) { face1 = 0; face2 = 3; }
      else if (id_node2 == pts[4]) { face1 = 3; face2 = 2; }
    } else if (id_node1 == pts[3]) {
      if      (id_node2 == pts[0]) { face1 = 0; face2 = 4; }
      else if (id_node2 == pts[2]) { face1 = 3; face2 = 0; }
      else if (id_node2 == pts[4]) { face1 = 4; face2 = 3; }
    } else if (id_node1 == pts[4]) {
      if      (id_node2 == pts[0]) { face1 = 4; face2 = 1; }
      else if (id_node2 == pts[1]) { face1 = 1; face2 = 2; }
      else if (id_node2 == pts[2]) { face1 = 2; face2 = 3; }
      else if (id_node2 == pts[3]) { face1 = 3; face2 = 4; }
    }
  } else {
    EG_BUG;
  }

  if (face1 == -1 || face2 == -1) {
    EG_BUG;
  }
}

void PolyMesh::getSortedEdgeCells(vtkIdType id_node1, vtkIdType id_node2, QList<vtkIdType> &cells, bool &is_loop)
{
  cells.clear();
  vtkIdType id_start = -1;
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, i);
    if (isVolume(id_cell, m_Grid)) {
      if (m_Cell2PCell[id_cell] == -1) {
        vtkIdType num_pts, *pts;
        m_Grid->GetCellPoints(id_cell, num_pts, pts);
        for (int j = 0; j < num_pts; ++j) {
          if (pts[j] == id_node2) {
            id_start = id_cell;
            break;
          }
        }
      }
    }
  }
  if (id_start == -1) {
    EG_BUG;
  }
  vtkIdType id_cell_old = id_start;
  vtkIdType id_cell_new = id_start;
  cells.append(id_cell_new);
  bool finished = false;
  bool prepend = false;
  while (!finished) {
    do {
      int f1, f2;
      getFacesOfEdgeInsideCell(id_cell_new, id_node1, id_node2, f1, f2);
      vtkIdType id_neigh1 = m_Part.c2cGG(id_cell_new, f1);
      vtkIdType id_neigh2 = m_Part.c2cGG(id_cell_new, f2);
      if (id_neigh1 == id_cell_old) {
        id_cell_old = id_cell_new;
        id_cell_new = id_neigh2;
      } else {
        id_cell_old = id_cell_new;
        id_cell_new = id_neigh1;
      }
      if (prepend) {
        cells.prepend(id_cell_new);
      } else {
        cells.append(id_cell_new);
      }
    } while (id_cell_new != id_start  &&  !isSurface(id_cell_new, m_Grid)  &&  m_Cell2PCell[id_cell_new] == -1);
    if ((isSurface(id_cell_new, m_Grid) || m_Cell2PCell[id_cell_new] != -1) && !prepend) {
      id_cell_new = cells[0];
      id_cell_old = cells[1];
      prepend = true;
    } else {
      finished = true;
    }
  }

  // check orientation (has to be from id_node1 to id_node2
  vec3_t x1, x2;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  vec3_t c = 0.5*(x1 + x2);
  vec3_t e = x2 - x1;
  QList<vtkIdType>::iterator iter_cell = cells.begin();
  vtkIdType id_cell1 = *iter_cell;
  ++iter_cell;
  bool reorientate = false;
  while (iter_cell != cells.end()) {
    vtkIdType id_cell2 = *iter_cell;
    vec3_t xc1 = cellCentre(m_Grid, id_cell1);
    vec3_t xc2 = cellCentre(m_Grid, id_cell2);
    vec3_t u = xc1 - c;
    vec3_t v = xc2 - c;
    vec3_t n = u.cross(v);
    if (n*e > 0 && reorientate) {
      EG_BUG;
    }
    if (n*e <= 0) {
      reorientate = true;
    }
    id_cell1 = id_cell2;
    ++iter_cell;
  }
  if (reorientate) {
    QList<vtkIdType> cells_old = cells;
    cells.clear();
    foreach (vtkIdType id_cell, cells_old) {
      cells.prepend(id_cell);
    }
  }

  // remove last cell for loops
  if (cells.size() > 1 && cells.first() == cells.last()) {
    cells.pop_back();
    is_loop = true;
  } else {
    is_loop = false;
  }
  /*
  if (prepend && (!isSurface(cells.first(), m_Grid) || !isSurface(cells.last(), m_Grid))) {
    EG_BUG;
  }
  if (!prepend && (isSurface(cells.first(), m_Grid) || isSurface(cells.last(), m_Grid))) {
    EG_BUG;
  }
  */
}

bool PolyMesh::isDualFace(vtkIdType id_face)
{
  vtkIdType id_cell = m_Part.getVolumeCell(id_face);
  if (m_Cell2PCell[id_cell] == -1) {
    return true;
  }
  return false;
}

void PolyMesh::getSortedPointFaces(vtkIdType id_node, int bc, QList<vtkIdType> &faces, bool &is_loop)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  faces.clear();
  vtkIdType id_start = -1;
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node, i);
    if (isSurface(id_cell, m_Grid)) {
      if (cell_code->GetValue(id_cell) == bc && isDualFace(id_cell)) {
        id_start = id_cell;
        break;
      }
    }
  }
  if (id_start == -1) {
    return;
  }
  vtkIdType id_face_old = id_start;
  vtkIdType id_face_new = id_start;
  faces.append(id_face_new);
  bool finished = false;
  bool prepend = false;
  while (!finished) {
    do {
      QList<vtkIdType> id_neigh;
      for (int i = 0; i < m_Part.c2cGSize(id_face_new); ++i) {
        if (cellContainsNode(m_Grid, m_Part.c2cGG(id_face_new, i), id_node)) {
          id_neigh.append(m_Part.c2cGG(id_face_new, i));
        }
      }
      if (id_neigh[0] == id_face_old) {
        id_face_old = id_face_new;
        id_face_new = id_neigh[1];
      } else {
        id_face_old = id_face_new;
        id_face_new = id_neigh[0];
      }
      if (cell_code->GetValue(id_face_new) == bc  &&  isDualFace(id_face_new)) {
      //if (cell_code->GetValue(id_face_new) == bc) {
        if (prepend) {
          faces.prepend(id_face_new);
        } else {
          faces.append(id_face_new);
        }
      }
    } while (id_face_new != id_start  &&  cell_code->GetValue(id_face_new) == bc  &&  isDualFace(id_face_new));
    if ((cell_code->GetValue(id_face_new) != bc || !isDualFace(id_face_new))  &&  !prepend) {
      id_face_old = id_face_new;
      id_face_new = faces[0];
      if (faces.size() > 1) {
        id_face_old = faces[1];
      }
      prepend = true;
    } else {
      finished = true;
    }
  }

  // remove last face for loops
  if (faces.size() > 1 && faces.first() == faces.last()) {
    faces.pop_back();
    is_loop = true;
  } else {
    is_loop = false;
  }

}

void PolyMesh::findPolyCells()
{
  m_Cell2PCell.fill(-1, m_Grid->GetNumberOfCells());
  m_Node2PCell.fill(-1, m_Grid->GetNumberOfPoints());
  m_NumPolyCells = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (!isHexCoreNode(id_node)) {
      bool tetra_or_pyramid = false;
      for (int j = 0; j < m_Part.n2cGSize(id_node); ++j) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, j);
        vtkIdType type_cell = m_Grid->GetCellType(id_cell);
        if (type_cell == VTK_TETRA || type_cell == VTK_PYRAMID) {
          tetra_or_pyramid = true;
        }
      }
      if (tetra_or_pyramid) {
        m_Node2PCell[id_node] = m_NumPolyCells;
        ++m_NumPolyCells;
      }
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if (type_cell == VTK_WEDGE || type_cell == VTK_HEXAHEDRON) {
      m_Cell2PCell[id_cell] = m_NumPolyCells;
      ++m_NumPolyCells;
    }
  }
}

void PolyMesh::createFace(QList<node_t> nodes, int owner, int neighbour, vec3_t ref_vec, int bc)
{
  if (owner > neighbour && neighbour != -1) {
    swap(owner, neighbour);
    ref_vec *= -1;
  }
  face_t face(nodes.size(), owner, neighbour, ref_vec, bc);
  for (int i = 0; i < nodes.size(); ++i) {
    int idx = m_Nodes.insert(nodes[i]);
    face.node[i] = idx;
  }
  m_Faces.append(face);
}

void PolyMesh::createCornerFace(vtkIdType id_cell, int i_face, vtkIdType id_node)
{
  QList<vtkIdType> edge_nodes;
  QVector<vtkIdType> face_nodes;
  getFaceOfCell(m_Grid, id_cell, i_face, face_nodes);
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell, num_pts, pts);
  //vtkIdType id_node = pts[i_node];
  for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
    vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
    if (face_nodes.contains(id_neigh)) {
      edge_nodes.append(id_neigh);
    }
  }
  if (edge_nodes.size() != 2) {
    EG_BUG;
  }
  int owner = m_Cell2PCell[id_cell];
  if (owner == -1) {
    EG_BUG;
  }
  int neighbour = m_Node2PCell[id_node];
  int bc = 0;
  if (neighbour == -1) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    vtkIdType id_face = m_Part.c2cGG(id_cell, i_face);
    if (id_face == -1) {
      EG_BUG;
    }
    if (!isSurface(id_cell, m_Grid)) {
      EG_BUG;
    }
    bc = cell_code->GetValue(id_face);
  }
  QList<node_t> nodes;
  nodes.append(node_t(id_node));
  nodes.append(node_t(id_node, edge_nodes[0]));
  nodes.append(node_t(face_nodes));
  nodes.append(node_t(id_node, edge_nodes[1]));
  vec3_t n = getNormalOfCell(m_Grid, id_cell, i_face);
  n.normalise();
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::createEdgeFace(vtkIdType id_node1, vtkIdType id_node2)
{
  // check id additional edge node needs to be created

  bool dual1 = false;
  bool dual2 = false;
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    if (m_Cell2PCell[m_Part.n2cGG(id_node1,i)] != -1) {
      dual1 = true;
      break;
    }
  }
  for (int i = 0; i < m_Part.n2cGSize(id_node2); ++i) {
    if (m_Cell2PCell[m_Part.n2cGG(id_node2,i)] != -1) {
      dual2 = true;
      break;
    }
  }
  bool add_edge_node = dual1 && dual2;

  /*
  if (id_node1 == 109 && id_node2 == 119) {
    cout << "break";
  }
  */

  QList<vtkIdType> cells;
  bool loop = false;
  getSortedEdgeCells(id_node1, id_node2, cells, loop);
  int owner     = m_Node2PCell[id_node1];
  int neighbour = m_Node2PCell[id_node2];
  if (owner == -1) {
    EG_BUG;
  }
  if (neighbour == -1) {
    EG_BUG;
  }
  if (owner > neighbour) {
    swap(id_node1, id_node2);
    swap(owner, neighbour);
  }
  QList<node_t> nodes;
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    node_t node;
    vtkIdType num_pts1, *pts1;
    m_Grid->GetCellPoints(cells[i_cells], num_pts1, pts1);
    if (m_Cell2PCell[cells[i_cells]] != -1) {
      int i_cells2 = 0;
      if (i_cells == 0) {
        i_cells2 = 1;
      } else if (i_cells == cells.size() - 1) {
        i_cells2 = cells.size() - 2;
      } else {
        EG_BUG;
      }
      vtkIdType num_pts2, *pts2;
      m_Grid->GetCellPoints(cells[i_cells2], num_pts2, pts2);
      QSet<vtkIdType> p1, p2;
      for (int i_pts1 = 0; i_pts1 < num_pts1; ++i_pts1) {
        p1.insert(pts1[i_pts1]);
      }
      for (int i_pts2 = 0; i_pts2 < num_pts2; ++i_pts2) {
        p2.insert(pts2[i_pts2]);
      }
      QSet<vtkIdType> face_nodes = p1.intersect(p2);
      node.id.resize(face_nodes.size());
      int i = 0;
      foreach (vtkIdType id_node, face_nodes) {
        node.id[i] = id_node;
        ++i;
      }
    } else {
      node.id.resize(num_pts1);
      for (int i_pts1 = 0; i_pts1 < num_pts1; ++i_pts1) {
        node.id[i_pts1] = pts1[i_pts1];
      }
    }
    qSort(node.id);
    nodes.append(node);
  }
  if (!loop) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    if (cell_code->GetValue(cells.first()) != cell_code->GetValue(cells.last()) || add_edge_node) {
      nodes.append(node_t(id_node1, id_node2));
    }
  }
  vec3_t x1, x2;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  createFace(nodes, owner, neighbour, x2 - x1, 0);
}

void PolyMesh::createFaceFace(vtkIdType id_cell, int i_face)
{
  if (m_Cell2PCell[id_cell] == -1) {
    EG_BUG;
  }
  int owner = m_Cell2PCell[id_cell];
  int neighbour = -1;
  vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
  int bc = 0;
  if (id_neigh != -1) {
    if (isVolume(id_neigh, m_Grid)) {
      if (m_Cell2PCell[id_neigh] == -1) {
        EG_BUG;
      }
      neighbour = m_Cell2PCell[id_neigh];
    } else {
      EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
      bc = cell_code->GetValue(id_neigh);
    }
  }
  QVector<vtkIdType> tmp_node_ids;
  getFaceOfCell(m_Grid, id_cell, i_face, tmp_node_ids);
  vec3_t n = getNormalOfCell(m_Grid, id_cell, i_face);
  n.normalise();
  QVector<vtkIdType> node_ids(tmp_node_ids.size() + 1);
  for (int i = 0; i < tmp_node_ids.size(); ++i) {
    node_ids[i] = tmp_node_ids[i];
  }
  node_ids[node_ids.size() - 1] = node_ids[0];
  QList<node_t> nodes;
  for (int i = 0; i < node_ids.size() - 1; ++i) {
    nodes.append(node_t(node_ids[i]));
    if (m_Node2PCell[node_ids[i]] != -1 && m_Node2PCell[node_ids[i+1]] != -1) {
      nodes.append(node_t(node_ids[i], node_ids[i+1]));
    }
  }
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::createPointFace(vtkIdType id_node, int bc)
{
  bool is_loop;
  QList<vtkIdType> faces;
  getSortedPointFaces(id_node, bc, faces, is_loop);
  if (faces.size() == 0) {
    return;
  }
  QList<node_t> nodes;
  vec3_t n(0,0,0);
  if (faces.size() == 0) {
    EG_BUG;
  }
  foreach (vtkIdType id_face, faces) {
    node_t node;
    vtkIdType num_pts, *pts;
    m_Grid->GetCellPoints(id_face, num_pts, pts);
    node.id.resize(num_pts);
    for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
      node.id[i_pts] = pts[i_pts];
    }
    qSort(node.id);
    nodes.append(node);
    n += GeometryTools::cellNormal(m_Grid, id_face);
  }
  n.normalise();
  if (!is_loop) {
    bool prepend = false;
    vtkIdType id_face = faces.last();
    while (id_face != -1) {
      bool found = false;
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_face, num_pts, pts);
      QList<vtkIdType> id_neigh_node;
      vtkIdType id_neigh = -1;
      for (int i = 0; i < m_Part.c2cGSize(id_face); ++i) {
        id_neigh = m_Part.c2cGG(id_face, i);
        if (id_neigh == -1) {
          EG_BUG;
        }
        if (!isSurface(id_neigh, m_Grid)) {
          EG_BUG;
        }
        if (cellContainsNode(m_Grid, id_neigh, id_node)) {
          if (!faces.contains(id_neigh)) {
            if (found) {
              EG_BUG;
            }
            for (int j = 0; j < num_pts; ++j) {
              if (pts[j] != id_node && cellContainsNode(m_Grid, id_neigh, pts[j])) {
                id_neigh_node.append(pts[j]);
              }
            }
          }
        }
      }

      if (id_neigh_node.size() == 0) {
        EG_BUG;
      }

      if (id_neigh_node.size() == 1) {
        if (prepend) {
          nodes.prepend(node_t(id_node, id_neigh_node[0]));
          id_face = -1;
        } else {
          nodes.append(node_t(id_node, id_neigh_node[0]));
          prepend = true;
          id_face = faces.first();
        }
      } else {
        if (id_neigh_node.size() > 2) {
          EG_BUG;
        }
        if (faces.size() != 1) {
          EG_BUG;
        }
        nodes.prepend(node_t(id_node, id_neigh_node[0]));
        /*
        node_t node;
        node.id.resize(num_pts);
        for (int i = 0; i < num_pts; ++i) {
          node.id[i] = pts[i];
        }
        qSort(node.id);
        nodes.append(node);
        */
        nodes.append(node_t(id_node, id_neigh_node[1]));
        id_face = -1;
      }
    }
    nodes.prepend(node_t(id_node));
  }
  int owner     = m_Node2PCell[id_node];
  int neighbour = -1;
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::createNodesAndFaces()
{
  m_Nodes.resize(m_Grid->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  // create all prismatic elements (hexes and prisms)
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Cell2PCell[id_cell] != -1) {
      for (int i_face = 0; i_face < m_Part.c2cGSize(id_cell); ++i_face) {
        vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
        if (id_neigh == -1) {
          EG_BUG;
        }
        bool create_corner_faces = false;
        if (m_Cell2PCell[id_neigh] == -1 && !isSurface(id_neigh, m_Grid)) {
          create_corner_faces = true;
        }
        if (create_corner_faces) {
          QVector<vtkIdType> face_nodes;
          getFaceOfCell(m_Grid, id_cell, i_face, face_nodes);
          foreach (vtkIdType id_node, face_nodes) {
            if (m_Node2PCell[id_node] == -1) {
              EG_BUG;
            }
            createCornerFace(id_cell, i_face, id_node);
          }
        } else {
          if (id_neigh > id_cell || isSurface(id_neigh, m_Grid)) {
            createFaceFace(id_cell, i_face);
          }
        }
      }
    }
  }

  // create all dual cells
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_Node2PCell[id_node] != -1) {
      for (int i_neigh = 0; i_neigh < m_Part.n2nGSize(id_node); ++i_neigh) {
        vtkIdType id_neigh = m_Part.n2nGG(id_node, i_neigh);
        if (m_Node2PCell[id_neigh] != -1 && id_neigh > id_node) {
          createEdgeFace(id_node, id_neigh);
        }
      }
      QSet<int> bcs;
      for (int i_cell = 0; i_cell < m_Part.n2cGSize(id_node); ++i_cell) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i_cell);
        if (isSurface(id_cell, m_Grid)) {
          bcs.insert(cell_code->GetValue(id_cell));
        }
      }
      foreach (int bc, bcs) {
        createPointFace(id_node, bc);
      }
    }
  }

  // compute node coordinates
  QVector<node_t> nodes;
  m_Nodes.getQVector(nodes);
  m_Points.resize(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    if (nodes[i].id.size() == 0) {
      EG_BUG;
    }
    m_Points[i] = vec3_t(0,0,0);
    //cout << i << ": ";
    foreach (vtkIdType id, nodes[i].id) {
      vec3_t x;
      //cout << id << " ";
      m_Grid->GetPoint(id, x.data());
      m_Points[i] += x;
    }
    //cout << endl;
    m_Points[i] *= 1.0/nodes[i].id.size();
  }

  qSort(m_Faces);

  QSet<int> bcs;
  foreach (face_t face, m_Faces) {
    if (face.bc != 0) {
      bcs.insert(face.bc);
    }
  }
  m_BCs.resize(bcs.size());
  qCopy(bcs.begin(), bcs.end(), m_BCs.begin());

}

vec3_t PolyMesh::faceNormal(int i)
{
  int N = m_Faces[i].node.size();
  QVector<vec3_t> x(N + 1);
  vec3_t xc(0,0,0);
  for (int j = 0; j < N; ++j) {
    x[j] = m_Points[m_Faces[i].node[j]];
    xc += x[j];
  }
  x[N] = x[0];
  xc *= 1.0/N;
  vec3_t n(0,0,0);
  for (int j = 0; j < N; ++j) {
    vec3_t u = x[j] - xc;
    vec3_t v = x[j+1] - xc;
    n += 0.5*u.cross(v);
  }
  return n;
}

void PolyMesh::checkFaceOrientation()
{
  for (int i = 0; i < m_Faces.size(); ++i) {
    int N = m_Faces[i].node.size();
    vec3_t n = faceNormal(i);
    n.normalise();
    if (n*m_Faces[i].ref_vec < 0) {
      QVector<int> old_node = m_Faces[i].node;
      for (int j = 0; j < N; ++j) {
        m_Faces[i].node[j] = old_node[N-1-j];
      }
    }
  }
}

