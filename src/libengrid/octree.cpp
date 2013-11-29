// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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

#include "octree.h"
#include "geometrytools.h"
#include "deletestraynodes.h"
#include "guimainwindow.h"

OctreeCell::OctreeCell()
{
  for (int i = 0; i < 8; ++i) {
    m_Node[i] = -1;
    m_Child[i] = -1;
  }
  for (int i = 0; i < 6; ++i) {
    m_Neighbour[i] = -1;
  }
  m_Parent = -1;
  m_Level = 0;
}

int OctreeCell::findNode(int i_node)
{
  int i_local = -1;
  for (int i = 0; i < 8; ++i) {
    if (m_Node[i] == i_node) {
      i_local = i;
      break;
    }
  }
  return i_local;
}

int OctreeCell::getEdgeNode(Octree* octree, int n1, int n2, bool this_cell_only)
{
  vec3_t x1 = octree->getNodePosition(m_Node[n1]);
  vec3_t x2 = octree->getNodePosition(m_Node[n2]);
  vec3_t xm = 0.5*(x1 + x2);
  int edge_node = -1;
  if (hasChildren()) {
    if (n1 > n2) {
      swap(n1, n2);
    }
    if (n1 == 0) {
      if (n2 == 1 || n2 == 2 || n2 == 4) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 1) {
      if (n2 == 3 || n2 == 5) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 2) {
      if (n2 == 3 || n2 == 6) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 3) {
      if (n2 == 7) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 4) {
      if (n2 == 5 || n2 == 6) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 5) {
      if (n2 == 7) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    } else if (n1 == 6) {
      if (n2 == 7) {
        edge_node = octree->m_Cells[m_Child[n1]].m_Node[n2];
      }
    }
  } else if (!this_cell_only) {
    QSet<int> cells1 = QSet<int>::fromList(octree->m_Node2Cell[m_Node[n1]]);
    QSet<int> cells2 = QSet<int>::fromList(octree->m_Node2Cell[m_Node[n2]]);
    QSet<int> cells = cells1.intersect(cells2);
    foreach (int i_cell, cells) {
      if (octree->m_Cells[i_cell].m_Level == m_Level) {
        int nn1 = octree->m_Cells[i_cell].findNode(m_Node[n1]);
        int nn2 = octree->m_Cells[i_cell].findNode(m_Node[n2]);
        edge_node = octree->m_Cells[i_cell].getEdgeNode(octree, nn1, nn2, true);
        if (edge_node != -1) {
          break;
        }
      }
    }
  }
  return edge_node;
}

void OctreeCell::getFaceNodes(int i, Octree* octree, QVector<int>& face_nodes, bool reverse)
{
  //face_nodes.resize(4);
  QList<int> nodes;
  int edge_node;
  if (i == 0) {
    nodes.push_back(m_Node[0]);
    edge_node = getEdgeNode(octree, 0,4); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[4]);
    edge_node = getEdgeNode(octree, 4,6); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[6]);
    edge_node = getEdgeNode(octree, 6,2); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[2]);
    edge_node = getEdgeNode(octree, 2,0); if (edge_node != -1) nodes.push_back(edge_node);
  } else if (i == 1) {
    nodes.push_back(m_Node[1]);
    edge_node = getEdgeNode(octree, 1,3); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[3]);
    edge_node = getEdgeNode(octree, 3,7); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[7]);
    edge_node = getEdgeNode(octree, 7,5); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[5]);
    edge_node = getEdgeNode(octree, 5,1); if (edge_node != -1) nodes.push_back(edge_node);
  } else if (i == 2) {
    nodes.push_back(m_Node[0]);
    edge_node = getEdgeNode(octree, 0,1); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[1]);
    edge_node = getEdgeNode(octree, 1,5); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[5]);
    edge_node = getEdgeNode(octree, 5,4); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[4]);
    edge_node = getEdgeNode(octree, 4,0); if (edge_node != -1) nodes.push_back(edge_node);
  } else if (i == 3) {
    nodes.push_back(m_Node[3]);
    edge_node = getEdgeNode(octree, 3,2); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[2]);
    edge_node = getEdgeNode(octree, 2,6); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[6]);
    edge_node = getEdgeNode(octree, 6,7); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[7]);
    edge_node = getEdgeNode(octree, 7,3); if (edge_node != -1) nodes.push_back(edge_node);
  } else if (i == 4) {
    nodes.push_back(m_Node[0]);
    edge_node = getEdgeNode(octree, 0,2); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[2]);
    edge_node = getEdgeNode(octree, 2,3); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[3]);
    edge_node = getEdgeNode(octree, 3,1); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[1]);
    edge_node = getEdgeNode(octree, 1,0); if (edge_node != -1) nodes.push_back(edge_node);
  } else if (i == 5) {
    nodes.push_back(m_Node[4]);
    edge_node = getEdgeNode(octree, 4,5); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[5]);
    edge_node = getEdgeNode(octree, 5,7); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[7]);
    edge_node = getEdgeNode(octree, 7,6); if (edge_node != -1) nodes.push_back(edge_node);
    nodes.push_back(m_Node[6]);
    edge_node = getEdgeNode(octree, 6,4); if (edge_node != -1) nodes.push_back(edge_node);
  }

  face_nodes.resize(nodes.size());
  if (reverse) {
    int i = face_nodes.size() - 1;
    foreach(int node, nodes) {
      face_nodes[i] = node;
      --i;
    }
  } else {
    int i = 0;
    foreach(int node, nodes) {
      face_nodes[i] = node;
      ++i;
    }
  }
}

void OctreeCell::getFaceNodes(int i, Octree* octree, QVector<QVector<int> >& face_nodes, bool reverse)
{
  if (hasChildren()) {
    face_nodes.resize(4);
    if (i == 0) {
      octree->m_Cells[m_Child[0]].getFaceNodes(0, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[2]].getFaceNodes(0, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[4]].getFaceNodes(0, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[6]].getFaceNodes(0, octree, face_nodes[3], reverse);
    } else if (i == 1) {
      octree->m_Cells[m_Child[1]].getFaceNodes(1, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[3]].getFaceNodes(1, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[5]].getFaceNodes(1, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[7]].getFaceNodes(1, octree, face_nodes[3], reverse);
    } else if (i == 2) {
      octree->m_Cells[m_Child[0]].getFaceNodes(2, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[1]].getFaceNodes(2, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[4]].getFaceNodes(2, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[5]].getFaceNodes(2, octree, face_nodes[3], reverse);
    } else if (i == 3) {
      octree->m_Cells[m_Child[2]].getFaceNodes(3, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[3]].getFaceNodes(3, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[6]].getFaceNodes(3, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[7]].getFaceNodes(3, octree, face_nodes[3], reverse);
    } else if (i == 4) {
      octree->m_Cells[m_Child[0]].getFaceNodes(4, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[1]].getFaceNodes(4, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[2]].getFaceNodes(4, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[3]].getFaceNodes(4, octree, face_nodes[3], reverse);
    } else if (i == 5) {
      octree->m_Cells[m_Child[4]].getFaceNodes(5, octree, face_nodes[0], reverse);
      octree->m_Cells[m_Child[5]].getFaceNodes(5, octree, face_nodes[1], reverse);
      octree->m_Cells[m_Child[6]].getFaceNodes(5, octree, face_nodes[2], reverse);
      octree->m_Cells[m_Child[7]].getFaceNodes(5, octree, face_nodes[3], reverse);
    }
  } else {
    face_nodes.resize(1);
    getFaceNodes(i, octree, face_nodes[0], reverse);
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Octree::Octree()
{
  m_Origin = vec3_t(0,0,0);
  setBase(vec3_t(1,0,0), vec3_t(0,1,0), vec3_t(0,0,1));
  m_Nodes.resize(8);
  m_Cells.resize(1);
  m_ToRefine.fill(false, 1);
  setBounds(vec3_t(0,0,0), vec3_t(1,1,1));
  setSmoothTransitionOn();
  setMaxCells(1000000);
}

void Octree::markToRefine(int cell)
{
  if (!hasChildren(cell)) {
    m_ToRefine[cell] = true;
  }
}

void Octree::markAllToRefine()
{
  for (int cell = 0; cell < m_Cells.size(); ++cell) {
    markToRefine(cell);
  }
}

void Octree::mergeNodes_identifyDuplicates()
{
  m_SameNodes.resize(m_Nodes.size());
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    m_SameNodes[i_nodes] = i_nodes;
  }


  int max_level = 0;
  for (int i_cell = 0; i_cell < m_Cells.size(); ++i_cell) {
    max_level = max(max_level, m_Cells[i_cell].m_Level);
  }
  buildNode2Cell();

  for (int i_node = 0; i_node < getNumNodes(); ++i_node) {
    QSet<int> canditates;
    double tol = 1e99;
    foreach (int i_cell, m_Node2Cell[i_node]) {
      tol = min(tol, 0.01*min(getDx(i_cell), min(getDy(i_cell), getDz(i_cell))));
      for (int i = 0; i < 8; ++i) {
        int j_node = m_Cells[i_cell].m_Node[i];
        foreach (int j_cell, m_Node2Cell[j_node]) {
          for (int j = 0; j < 8; ++j) {
            canditates.insert(m_Cells[j_cell].m_Node[j]);
          }
        }
      }
    }
    if (tol > 1e98) {
      EG_BUG;
    }
    foreach (int j_node, canditates) {
      if (i_node != j_node) {
        double distance = (m_Nodes[i_node].m_Position - m_Nodes[j_node].m_Position).abs();
        if (distance < tol) {
          if (i_node > j_node) {
            m_SameNodes[i_node] = j_node;
          } else {
            m_SameNodes[j_node] = i_node;
          }
        }
      }
    }
  }

  QVector<bool> is_dup(m_Nodes.size(), false);
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {

    int i_same_node = i_nodes;
    while (i_same_node != m_SameNodes[i_same_node]) {
      i_same_node = m_SameNodes[i_same_node];
    }
    m_SameNodes[i_nodes] = i_same_node;
    if (m_SameNodes[i_nodes] != i_nodes) {
      is_dup[m_SameNodes[i_nodes]] = true;
    }
    if (m_SameNodes[m_SameNodes[i_nodes]] != m_SameNodes[i_nodes]) {
      EG_BUG;
    }
  }
  int N = 0;
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    if (is_dup[i_nodes]) {
      ++N;
    }
  }
}

void Octree::mergeNodes_compactNodes()
{
  QVector<int> offset(m_Nodes.size());
  int last_offset = 0;
  double max_dist = 0;
  for (int i = 0; i < m_Nodes.size(); ++i) {
    if (m_SameNodes[i] != i) {
      max_dist = max(max_dist, (m_Nodes[i].m_Position - m_Nodes[m_SameNodes[i]].m_Position).abs());
      if (m_SameNodes[i] > i) {
        EG_BUG;
      }
      ++last_offset;
    }
    offset[i] = last_offset;
  }

  for (int i = 0; i < m_Nodes.size(); ++i) {
    if (m_SameNodes[i] != i) {
      m_SameNodes[i] -= offset[m_SameNodes[i]];
    } else {
      m_SameNodes[i] -= offset[i];
    }
  }

  for (int i = 0; i < m_Nodes.size(); ++i) {
    if (m_SameNodes[i] > i) {
      EG_BUG;
    }
  }

  QVector<int> copy_nodes(m_Nodes.size() - last_offset);
  for (int i = 0; i < m_Nodes.size(); ++i) {
    copy_nodes[m_SameNodes[i]] = i;
  }

  for (int i = 0; i < m_Nodes.size() - last_offset; ++i) {
    m_Nodes[i] = m_Nodes[copy_nodes[i]];
  }

  m_Nodes.remove(m_Nodes.size() - last_offset, last_offset);

}

void Octree::mergeNodes_updateCells()
{
  for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
    for (int j = 0; j < 8; ++j) {
      m_Cells[i_cells].m_Node[j] = m_SameNodes[m_Cells[i_cells].m_Node[j]];
    }
  }
}

void Octree::checkNeighbours()
{
  int Nerr = 0;
  int Nno = 0;
  for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
    for (int i_faces = 0; i_faces < 6; ++i_faces) {
      int i_neigh_cells = m_Cells[i_cells].m_Neighbour[i_faces];
      if (i_neigh_cells != -1) {
        int i_back_faces = i_faces;
        if (i_faces % 2 == 0) {
          ++i_back_faces;
        } else {
          --i_back_faces;
        }
        if (m_Cells[i_cells].m_Level == m_Cells[i_neigh_cells].m_Level) {
          if (m_Cells[i_neigh_cells].m_Neighbour[i_back_faces] != i_cells) {
            cout << "neighbour error: " << i_cells << ',' << i_neigh_cells << ',' << i_faces << ',' << i_back_faces << ','
                 << getCellCentre(i_cells) << ','  << getCellCentre(i_neigh_cells) << endl;
            ++Nerr;
          }
        }
      } else {
        cout << "no neighbour: " << i_cells << ',' << i_faces << ',' << getFaceCentre(i_cells, i_faces) << endl;
        ++Nno;
      }
    }
  }
  cout << Nerr << " errors and " << Nno << " faces without neighbour" << endl;
}

void Octree::buildNode2Cell()
{
  m_Node2Cell.clear();
  m_Node2Cell.fill(QList<int>(), getNumNodes());
  for (int i_cell = 0; i_cell < getNumCells(); ++i_cell) {
    for (int j = 0; j < 8; ++j) {
      m_Node2Cell[getNode(i_cell, j)].append(i_cell);
    }
  }
  for (int i_node = 0; i_node < getNumNodes(); ++i_node) {
    int max_level = 1;
    foreach (int i_cell, m_Node2Cell[i_node]) {
      max_level = max(max_level, m_Cells[i_cell].m_Level);
    }
  }
}

void Octree::mergeNodes()
{
  mergeNodes_identifyDuplicates();
  mergeNodes_compactNodes();
  mergeNodes_updateCells();
}

void Octree::setOrigin(vec3_t x0)
{
  m_Origin = x0;
}

void Octree::setBase(vec3_t g1, vec3_t g2, vec3_t g3)
{
  m_Base.column(0, g1);
  m_Base.column(1, g2);
  m_Base.column(2, g3);
  m_InvBase = m_Base.inverse();
}

void Octree::setBounds(vec3_t corner1, vec3_t corner2, int num_i, int num_j, int num_k)
{
  vec3_t g1(m_Base[0][0], m_Base[1][0], m_Base[2][0]);
  vec3_t g2(m_Base[0][1], m_Base[1][1], m_Base[2][1]);
  vec3_t g3(m_Base[0][2], m_Base[1][2], m_Base[2][2]);

  int num_cells = num_i*num_j*num_k;
  int num_nodes = (num_i + 1)*(num_j + 1)*(num_k + 1);
  m_Cells.resize(num_cells);
  m_Nodes.resize(num_nodes);
  m_Corner1 = corner1;
  m_Corner2 = corner2;
  m_Dx = fabs(corner1[0] - corner2[0])/num_i;
  m_Dy = fabs(corner1[1] - corner2[1])/num_j;
  m_Dz = fabs(corner1[2] - corner2[2])/num_k;

  for (int i = 0; i < num_i + 1; ++i) {
    for (int j = 0; j < num_j + 1; ++j) {
      for (int k = 0; k < num_k + 1; ++k) {
        int idx_node = i*(num_j + 1)*(num_k + 1) + j*(num_k + 1) + k;
        m_Nodes[idx_node].setPosition(m_Corner1 + i*m_Dx*g1 + j*m_Dy*g2 + k*m_Dz*g3);
      }
    }
  }
  for (int i = 0; i < num_i; ++i) {
    for (int j = 0; j < num_j; ++j) {
      for (int k = 0; k < num_k; ++k) {
        int idx_cell  = i*num_j*num_k + j*num_k + k;
        int idx_node0 = (i + 0)*(num_j + 1)*(num_k + 1) + (j + 0)*(num_k + 1) + (k + 0);
        int idx_node1 = (i + 1)*(num_j + 1)*(num_k + 1) + (j + 0)*(num_k + 1) + (k + 0);
        int idx_node2 = (i + 0)*(num_j + 1)*(num_k + 1) + (j + 1)*(num_k + 1) + (k + 0);
        int idx_node3 = (i + 1)*(num_j + 1)*(num_k + 1) + (j + 1)*(num_k + 1) + (k + 0);
        int idx_node4 = (i + 0)*(num_j + 1)*(num_k + 1) + (j + 0)*(num_k + 1) + (k + 1);
        int idx_node5 = (i + 1)*(num_j + 1)*(num_k + 1) + (j + 0)*(num_k + 1) + (k + 1);
        int idx_node6 = (i + 0)*(num_j + 1)*(num_k + 1) + (j + 1)*(num_k + 1) + (k + 1);
        int idx_node7 = (i + 1)*(num_j + 1)*(num_k + 1) + (j + 1)*(num_k + 1) + (k + 1);
        m_Cells[idx_cell].m_Node[0] = idx_node0;
        m_Cells[idx_cell].m_Node[1] = idx_node1;
        m_Cells[idx_cell].m_Node[2] = idx_node2;
        m_Cells[idx_cell].m_Node[3] = idx_node3;
        m_Cells[idx_cell].m_Node[4] = idx_node4;
        m_Cells[idx_cell].m_Node[5] = idx_node5;
        m_Cells[idx_cell].m_Node[6] = idx_node6;
        m_Cells[idx_cell].m_Node[7] = idx_node7;
      }
    }
  }
}

int Octree::refineAll()
{
  int N1;
  int N2;
  int count = 0;
  for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
    if (m_Cells[i_cells].m_Child[0] != -1) {
      m_ToRefine[i_cells] = false;
    }
  }
  buildNode2Cell();
  do {
    N1 = 0;
    N2 = 0;
    for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
      if (m_ToRefine[i_cells]) {
        ++N1;
        if (m_SmoothTransition) {
          for (int i = 0; i < 8; ++i) {
            foreach (int neigh, m_Node2Cell[m_Cells[i_cells].m_Node[i]]) {
              if ((m_Cells[neigh].m_Level < m_Cells[i_cells].m_Level) && !m_ToRefine[neigh] && !hasChildren(neigh)) {
                markToRefine(neigh);
                ++N2;
              }
            }
          }
        }
      }
    }
    ++count;
    if (count > 100) {
      EG_BUG;
    }
  } while (N2 > 0);
  N2 = m_Cells.size();
  if (N2 + 8*N1 > m_MaxCells) {
    QString num;
    QString msg = "maximal number of cells exceeded\n";
    num.setNum(N2 + 8*N1);
    msg += num += " requested and ";
    num.setNum(m_MaxCells);
    msg += num + " allowed";
    EG_ERR_RETURN(msg);
  }
  m_Cells.insert(N2, 8*N1, OctreeCell());
  int Nrefine = N1;
  int new_node = m_Nodes.size();
  m_Nodes.insert(m_Nodes.size(), N1*19, OctreeNode());
  N1 = N2;
  for (int i_cells = 0; i_cells < N1; ++i_cells) {
    if (m_ToRefine[i_cells]) {
      if (hasChildren(i_cells)) {
        EG_BUG;
      }
      int nn[8];
      nn[0] = m_Cells[i_cells].m_Node[0];
      nn[1] = m_Cells[i_cells].m_Node[1];
      nn[2] = m_Cells[i_cells].m_Node[2];
      nn[3] = m_Cells[i_cells].m_Node[3];
      nn[4] = m_Cells[i_cells].m_Node[4];
      nn[5] = m_Cells[i_cells].m_Node[5];
      nn[6] = m_Cells[i_cells].m_Node[6];
      nn[7] = m_Cells[i_cells].m_Node[7];
      // create new nodes
      int ne[12];
      ne[0] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position);
      ne[1] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position);
      ne[2] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position);
      ne[3] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position);
      ne[4] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position);
      ne[5] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      ne[6] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position);
      ne[7] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      ne[8] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position);
      ne[9] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position);
      ne[10] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position);
      ne[11] = new_node;
      m_Nodes[new_node++].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      int nf[6];
      nf[0] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[0]].m_Position + m_Nodes[ne[4]].m_Position
                                             + m_Nodes[ne[8]].m_Position  + m_Nodes[ne[10]].m_Position);
      nf[1] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[1]].m_Position + m_Nodes[ne[5]].m_Position
                                             + m_Nodes[ne[9]].m_Position  + m_Nodes[ne[11]].m_Position);
      nf[2] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[2]].m_Position + m_Nodes[ne[6]].m_Position
                                             + m_Nodes[ne[8]].m_Position  + m_Nodes[ne[9]].m_Position);
      nf[3] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[3]].m_Position + m_Nodes[ne[7]].m_Position
                                             + m_Nodes[ne[10]].m_Position + m_Nodes[ne[11]].m_Position);
      nf[4] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[0]].m_Position + m_Nodes[ne[1]].m_Position
                                             + m_Nodes[ne[2]].m_Position  + m_Nodes[ne[3]].m_Position);
      nf[5] = new_node;
      m_Nodes[new_node++].m_Position = 0.25*(  m_Nodes[ne[4]].m_Position + m_Nodes[ne[5]].m_Position
                                             + m_Nodes[ne[6]].m_Position  + m_Nodes[ne[7]].m_Position);
      int nv = new_node;
      m_Nodes[new_node++].m_Position = 1.0/6.0*(  m_Nodes[nf[0]].m_Position + m_Nodes[nf[1]].m_Position + m_Nodes[nf[2]].m_Position
                                                + m_Nodes[nf[3]].m_Position + m_Nodes[nf[4]].m_Position + m_Nodes[nf[5]].m_Position);

      for (int child = 0; child < 8; ++child) {
        m_Cells[i_cells].m_Child[child] = N2;
        m_Cells[N2].m_Parent = i_cells;
        m_Cells[N2].m_Level = m_Cells[i_cells].m_Level + 1;
        ++N2;
      }

      // child 0
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[0] = nn[0];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[1] = ne[2];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[2] = ne[0];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[3] = nf[4];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[4] = ne[8];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[5] = nf[2];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[6] = nf[0];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[7] = nv;
      // child 1
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[0] = ne[2];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[1] = nn[1];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[2] = nf[4];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[3] = ne[1];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[4] = nf[2];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[5] = ne[9];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[6] = nv;
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Node[7] = nf[1];
      // child 2
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[0] = ne[0];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[1] = nf[4];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[2] = nn[2];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[3] = ne[3];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[4] = nf[0];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[5] = nv;
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[6] = ne[10];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Node[7] = nf[3];
      // child 3
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[0] = nf[4];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[1] = ne[1];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[2] = ne[3];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[3] = nn[3];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[4] = nv;
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[5] = nf[1];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[6] = nf[3];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Node[7] = ne[11];
      // child 4
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[0] = ne[8];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[1] = nf[2];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[2] = nf[0];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[3] = nv;
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[4] = nn[4];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[5] = ne[6];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[6] = ne[4];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Node[7] = nf[5];
      // child 5
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[0] = nf[2];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[1] = ne[9];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[2] = nv;
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[3] = nf[1];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[4] = ne[6];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[5] = nn[5];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[6] = nf[5];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Node[7] = ne[5];
      // child 6
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[0] = nf[0];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[1] = nv;
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[2] = ne[10];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[3] = nf[3];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[4] = ne[4];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[5] = nf[5];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[6] = nn[6];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Node[7] = ne[7];
      // child 7
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[0] = nv;
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[1] = nf[1];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[2] = nf[3];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[3] = ne[11];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[4] = nf[5];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[5] = ne[5];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[6] = ne[7];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Node[7] = nn[7];

      // - - -
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[1] = m_Cells[i_cells].m_Child[1];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[3] = m_Cells[i_cells].m_Child[2];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[5] = m_Cells[i_cells].m_Child[4];
      // + - -
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[0] = m_Cells[i_cells].m_Child[0];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[3] = m_Cells[i_cells].m_Child[3];
      m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[5] = m_Cells[i_cells].m_Child[5];
      // - + -
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[1] = m_Cells[i_cells].m_Child[3];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[2] = m_Cells[i_cells].m_Child[0];
      m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[5] = m_Cells[i_cells].m_Child[6];
      // + + -
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[0] = m_Cells[i_cells].m_Child[2];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[2] = m_Cells[i_cells].m_Child[1];
      m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[5] = m_Cells[i_cells].m_Child[7];
      // - - +
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[1] = m_Cells[i_cells].m_Child[5];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[3] = m_Cells[i_cells].m_Child[6];
      m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[4] = m_Cells[i_cells].m_Child[0];
      // + - +
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[0] = m_Cells[i_cells].m_Child[4];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[3] = m_Cells[i_cells].m_Child[7];
      m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[4] = m_Cells[i_cells].m_Child[1];
      // - + +
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[1] = m_Cells[i_cells].m_Child[7];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[2] = m_Cells[i_cells].m_Child[4];
      m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[4] = m_Cells[i_cells].m_Child[2];
      // + + +
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[0] = m_Cells[i_cells].m_Child[6];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[2] = m_Cells[i_cells].m_Child[5];
      m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[4] = m_Cells[i_cells].m_Child[3];
    }
  }

  for (int i_cells = 0; i_cells < N1; ++i_cells) {
    if (m_Cells[i_cells].m_Child[0] >= N1) {
      {
        int neigh = m_Cells[i_cells].m_Neighbour[0];
        if (neigh != -1) {
          if ((m_Cells[neigh].m_Child[0] != -1) && (m_Cells[i_cells].m_Level == m_Cells[neigh].m_Level)) {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[0] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[0] = m_Cells[neigh].m_Child[3];
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[0] = m_Cells[neigh].m_Child[5];
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[0] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[0] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[i_cells].m_Neighbour[1];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[1] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[1] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[1] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[1] = m_Cells[neigh].m_Child[6];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[1] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[i_cells].m_Neighbour[2];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[2] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[2] = m_Cells[neigh].m_Child[3];
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[2] = m_Cells[neigh].m_Child[6];
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[2] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[2] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[i_cells].m_Neighbour[3];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[3] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[3] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[3] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[3] = m_Cells[neigh].m_Child[5];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[3] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[i_cells].m_Neighbour[4];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[4] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[4] = m_Cells[neigh].m_Child[5];
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[4] = m_Cells[neigh].m_Child[6];
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[4] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[0]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[1]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[2]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[3]].m_Neighbour[4] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[i_cells].m_Neighbour[5];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[5] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[5] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[5] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[5] = m_Cells[neigh].m_Child[3];
          } else {
            m_Cells[m_Cells[i_cells].m_Child[4]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[5]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[6]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[i_cells].m_Child[7]].m_Neighbour[5] = neigh;
          }
        }
      }
    }
  }

  m_ToRefine.fill(false, m_Cells.size());
  mergeNodes();
  buildNode2Cell();
  return Nrefine;
}

void Octree::toVtkGridHangingNodes(vtkUnstructuredGrid *grid, bool)
{
  int N = 0;
  for (int i = 0; i < m_Cells.size(); ++i) {
    if (m_Cells[i].m_Child[0] == -1) {
      ++N;
    }
  }
  allocateGrid(grid, N, m_Nodes.size());
  //allocateGrid(grid, N, m_Nodes.size(), create_fields);
  for (int id_node = 0; id_node < m_Nodes.size(); ++id_node) {
    grid->GetPoints()->SetPoint(id_node, m_Nodes[id_node].m_Position.data());
  }
  EG_VTKDCC(vtkIntArray, cell_index, grid, "cell_code");
  for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
    OctreeCell cell = m_Cells[i_cells];
    if (cell.m_Child[0] == -1) {
      vtkIdType Npts = 8;
      vtkIdType pts[8];
      pts[0] = cell.m_Node[0];
      pts[1] = cell.m_Node[1];
      pts[2] = cell.m_Node[3];
      pts[3] = cell.m_Node[2];
      pts[4] = cell.m_Node[4];
      pts[5] = cell.m_Node[5];
      pts[6] = cell.m_Node[7];
      pts[7] = cell.m_Node[6];
      for (int i = 0; i < 8; ++i) {
        if (pts[i] < 0 || pts[i] >= m_Nodes.size()) {
          EG_BUG;
        }
      }
      vtkIdType id_cell = grid->InsertNextCell(VTK_HEXAHEDRON, Npts, pts);
      cell_index->SetValue(id_cell, i_cells);
    }
  }
}

int Octree::opposingFace(int i)
{
  if (i == 0) return 1;
  if (i == 1) return 0;
  if (i == 2) return 3;
  if (i == 3) return 2;
  if (i == 4) return 5;
  if (i == 5) return 4;
  EG_BUG;
  return -1;
}

void Octree::toVtkGridPolyhedral(vtkUnstructuredGrid *grid, bool create_fields)
{
  if (!m_SmoothTransition) {
    EG_BUG;
  }
  buildNode2Cell();
  vtkIdType num_new_cells = 0;
  for (int i_cell = 0; i_cell < getNumCells(); ++i_cell) {
    OctreeCell cell = m_Cells[i_cell];
    if (!cell.hasChildren()) { // only use cells which do not have children
      ++num_new_cells;
    }
  }
  allocateGrid(grid, num_new_cells, m_Nodes.size(), create_fields);
  for (vtkIdType id_node = 0; id_node < m_Nodes.size(); ++id_node) {
    grid->GetPoints()->SetPoint(id_node, m_Nodes[id_node].getPosition().data());
  }

  for (int i_cell = 0; i_cell < getNumCells(); ++i_cell) {
    OctreeCell cell = m_Cells[i_cell];
    if (!cell.hasChildren()) { // only use cells which do not have children
      QList<QVector<int> > all_faces;
      QVector<QVector<int> > faces;
      for (int i = 0; i < 6; ++i) {
        bool use_neighbour_faces = false;
        if (cell.getNeighbour(i) != -1) {
          OctreeCell neigh = m_Cells[cell.getNeighbour(i)];
          if (neigh.m_Level == cell.m_Level) {
            if (neigh.hasChildren()) {
              use_neighbour_faces = true;
            }
          }
        }
        if (use_neighbour_faces) {
          m_Cells[cell.getNeighbour(i)].getFaceNodes(opposingFace(i), this, faces, true);
        } else {
          cell.getFaceNodes(i, this, faces, false);
        }
        foreach (QVector<int> face, faces) {
          all_faces.push_back(face);
        }
      }
      if (all_faces.size() < 6) {
        EG_BUG;
      }
      bool simple_hex_cell = true;
      if (all_faces.size() > 6) {
        simple_hex_cell = false;
      }

      foreach (QVector<int> face, all_faces) {
        if (face.size() > 4) {
          simple_hex_cell = false;
          break;
        }
      }

      if (simple_hex_cell) {
        vtkIdType Npts = 8;
        vtkIdType pts[8];
        pts[0] = cell.m_Node[0];
        pts[1] = cell.m_Node[1];
        pts[2] = cell.m_Node[3];
        pts[3] = cell.m_Node[2];
        pts[4] = cell.m_Node[4];
        pts[5] = cell.m_Node[5];
        pts[6] = cell.m_Node[7];
        pts[7] = cell.m_Node[6];
        for (int i = 0; i < 8; ++i) {
          if (pts[i] < 0 || pts[i] >= m_Nodes.size()) {
            EG_BUG;
          }
        }
        vtkIdType id_cell = grid->InsertNextCell(VTK_HEXAHEDRON, Npts, pts);
      } else {
        vtkIdType num_ids = 1;
        foreach (QVector<int> face, all_faces) {
          num_ids += 1 + face.size();
        }
        EG_VTKSP(vtkIdList, ptids);
        ptids->SetNumberOfIds(num_ids);
        vtkIdType i_id = 0;
        ptids->SetId(i_id++, all_faces.size());
        foreach (QVector<int> face, all_faces) {
          ptids->SetId(i_id++, face.size());
          foreach (int i, face) {
            ptids->SetId(i_id++, i);
          }
        }
        grid->InsertNextCell(VTK_POLYHEDRON, ptids);
      }
    }
  }
}

void Octree::toVtkGridConforming(vtkUnstructuredGrid* grid, bool create_fields)
{
  if (!m_SmoothTransition) {
    EG_BUG;
  }
  buildNode2Cell();
  int num_new_nodes = 0;
  int num_pyramids = 0;
  int num_hexes = 0;
  int num_tetras = 0;
  for (int i_cell = 0; i_cell < getNumCells(); ++i_cell) {
    OctreeCell cell = m_Cells[i_cell];
    if (!cell.hasChildren()) { // only use cells which do not have children
      QList<QVector<int> > all_faces;
      QVector<QVector<int> > faces;
      for (int i = 0; i < 6; ++i) {
        bool use_neighbour_faces = false;
        if (cell.getNeighbour(i) != -1) {
          OctreeCell neigh = m_Cells[cell.getNeighbour(i)];
          if (neigh.m_Level == cell.m_Level) {
            if (neigh.hasChildren()) {
              use_neighbour_faces = true;
            }
          }
        }
        if (use_neighbour_faces) {
          m_Cells[cell.getNeighbour(i)].getFaceNodes(opposingFace(i), this, faces);
        } else {
          cell.getFaceNodes(i, this, faces, true);
        }
        foreach (QVector<int> face, faces) {
          all_faces.push_back(face);
        }
      }
      if (all_faces.size() < 6) {
        EG_BUG;
      }
      bool simple_hex_cell = true;
      if (all_faces.size() > 6) {
        simple_hex_cell = false;
      }

      foreach (QVector<int> face, all_faces) {
        if (face.size() > 4) {
          simple_hex_cell = false;
          break;
        }
      }

      if (simple_hex_cell) {
        ++num_hexes;
      } else {
        ++num_new_nodes;
        foreach (QVector<int> face, all_faces) {
          if (face.size() > 4) {
            num_tetras += face.size();
            ++num_new_nodes;
          } else {
            ++num_pyramids;
          }
        }
      }
    }
  }

  allocateGrid(grid, num_hexes + num_pyramids + num_tetras, m_Nodes.size() + num_new_nodes, create_fields);
  vtkIdType id_node = 0;
  for (int i = 0; i < m_Nodes.size(); ++i) {
    vec3_t x = m_Nodes[i].getPosition().data();
    grid->GetPoints()->SetPoint(id_node, m_Nodes[i].getPosition().data());
    ++id_node;
  }

  QVector<QList<node_t> > face_nodes(m_Cells.size());

  for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
    OctreeCell cell = m_Cells[i_cells];
    if (!cell.hasChildren()) {
      QList<QVector<int> > all_faces;
      QVector<QVector<int> > faces;
      for (int i = 0; i < 6; ++i) {
        bool use_neighbour_faces = false;
        if (cell.getNeighbour(i) != -1) {
          OctreeCell neigh = m_Cells[cell.getNeighbour(i)];
          if (neigh.m_Level == cell.m_Level) {
            if (neigh.m_Child[0] != -1) {
              use_neighbour_faces = true;
            }
          }
        }
        if (use_neighbour_faces) {
          m_Cells[cell.getNeighbour(i)].getFaceNodes(opposingFace(i), this, faces);
        } else {
          cell.getFaceNodes(i, this, faces, true);
        }
        foreach (QVector<int> face, faces) {
          all_faces.push_back(face);
        }
      }
      bool simple_hex_cell = true;
      if (all_faces.size() > 6) {
        simple_hex_cell = false;
      } else {
        foreach (QVector<int> face, all_faces) {
          if (face.size() > 4) {
            simple_hex_cell = false;
            break;
          }
        }
      }
      if (simple_hex_cell) {
        vtkIdType pts[8];
        pts[0] = cell.m_Node[0];
        pts[1] = cell.m_Node[1];
        pts[2] = cell.m_Node[3];
        pts[3] = cell.m_Node[2];
        pts[4] = cell.m_Node[4];
        pts[5] = cell.m_Node[5];
        pts[6] = cell.m_Node[7];
        pts[7] = cell.m_Node[6];
        grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
      } else {
        grid->GetPoints()->SetPoint(id_node, getCellCentre(i_cells).data());
        int id_cell_centre = id_node;
        ++id_node;
        foreach (QVector<int> face, all_faces) {
          if (face.size() == 4) {
            vtkIdType pts[5];
            for (int i = 0; i < 4; ++i) {
              pts[i] = face[i];
            }
            pts[4] = id_cell_centre;
            grid->InsertNextCell(VTK_PYRAMID, 5, pts);
          } else {
            vec3_t xf(0,0,0);
            vec3_t xc(0,0,0);
            for (int i = 0; i < face.size(); ++i) {
              xf += m_Nodes[face[i]].getPosition();
            }
            xf *= 1.0/face.size();
            double Atot = 0;
            double f13 = 1.0/3.0;
            for (int i = 0; i < face.size(); ++i) {
              vec3_t x1, x2;
              grid->GetPoint(face[i], x1.data());
              if (i < face.size()-1) {
                grid->GetPoint(face[i+1], x2.data());
              } else {
                grid->GetPoint(face[0], x2.data());
              }
              double A = GeometryTools::triArea(xf, x1, x2);
              Atot += A;
              xc += f13*A*(x1 + x2 + xf);
            }
            xc *= 1.0/Atot;

            int id_face_centre = -1;
            double l_crit = 0.001*min(getDx(i_cells), min(getDy(i_cells), getDz(i_cells)));
            QList<int> neighbours;
            for (int i = 0; i < 6; ++i) {
              if (cell.getNeighbour(i) != -1) {
                neighbours.append(cell.getNeighbour(i));
              }
            }
            foreach (int i_neigh, neighbours) {
              foreach (node_t N, face_nodes[i_neigh]) {
                if ((N.x - xc).abs() < l_crit) {
                  id_face_centre = N.id;
                }
              }
            }
            if (id_face_centre == -1) {
              grid->GetPoints()->SetPoint(id_node, xc.data());
              id_face_centre = id_node;
              ++id_node;
              node_t N;
              N.x = xc;
              N.id = id_face_centre;
              face_nodes[i_cells].append(N);
            }
            for (int i = 0; i < face.size(); ++i) {
              vtkIdType pts[4];
              pts[0] = face[i];
              if (i < face.size()-1) {
                pts[1] = face[i+1];
              } else {
                pts[1] = face[0];
              }
              pts[2] = id_face_centre;
              pts[3] = id_cell_centre;
              grid->InsertNextCell(VTK_TETRA, 4, pts);
            }
          }
        }        
      }
    }
  }

  DeleteStrayNodes del_stray;
  del_stray.setGrid(grid);
  del_stray.setAllCells();
  del_stray();

}

vec3_t Octree::getCellCentre(int cell)
{
  vec3_t x(0,0,0);
  for (int i = 0; i < 8; ++i) {
    x += m_Nodes[m_Cells[cell].m_Node[i]].m_Position;
  }
  x *= 1.0/8.0;
  return x;
}

vec3_t Octree::getFaceCentre(int i_cells, int i_faces)
{
  vec3_t x(0,0,0);
  const OctreeCell& cell = m_Cells[i_cells];
  if (i_faces == 0) {
    x += m_Nodes[cell.m_Node[0]].m_Position;
    x += m_Nodes[cell.m_Node[2]].m_Position;
    x += m_Nodes[cell.m_Node[4]].m_Position;
    x += m_Nodes[cell.m_Node[6]].m_Position;
  } else if (i_faces == 1) {
    x += m_Nodes[cell.m_Node[1]].m_Position;
    x += m_Nodes[cell.m_Node[3]].m_Position;
    x += m_Nodes[cell.m_Node[5]].m_Position;
    x += m_Nodes[cell.m_Node[7]].m_Position;
  } else if (i_faces == 2) {
    x += m_Nodes[cell.m_Node[0]].m_Position;
    x += m_Nodes[cell.m_Node[1]].m_Position;
    x += m_Nodes[cell.m_Node[4]].m_Position;
    x += m_Nodes[cell.m_Node[5]].m_Position;
  } else if (i_faces == 3) {
    x += m_Nodes[cell.m_Node[2]].m_Position;
    x += m_Nodes[cell.m_Node[3]].m_Position;
    x += m_Nodes[cell.m_Node[6]].m_Position;
    x += m_Nodes[cell.m_Node[7]].m_Position;
  } else if (i_faces == 4) {
    x += m_Nodes[cell.m_Node[0]].m_Position;
    x += m_Nodes[cell.m_Node[1]].m_Position;
    x += m_Nodes[cell.m_Node[2]].m_Position;
    x += m_Nodes[cell.m_Node[3]].m_Position;
  } else if (i_faces == 5) {
    x += m_Nodes[cell.m_Node[4]].m_Position;
    x += m_Nodes[cell.m_Node[5]].m_Position;
    x += m_Nodes[cell.m_Node[6]].m_Position;
    x += m_Nodes[cell.m_Node[7]].m_Position;
  }
  x *= 0.25;
  return x;
}

int Octree::findCell(vec3_t x)
{
  for (int i = 0; i < 3; ++i) {
    if ((x[i] < m_Corner1[i]) || (x[i] > m_Corner2[i])) {
      EG_ERR_RETURN("node outside of octree mesh");
    }
  }
  int i_cells = 0;
  while (getLevel(i_cells) == 0) {
    vec3_t x1 = m_Nodes[m_Cells[i_cells].m_Node[0]].getPosition();
    vec3_t x2 = m_Nodes[m_Cells[i_cells].m_Node[7]].getPosition();
    if (x[0] >= x1[0] && x[1] >= x1[1] && x[2] >= x1[2]  &&  x[0] <= x2[0] && x[1] <= x2[1] && x[2] <= x2[2]) {
      break;
    }
    ++i_cells;
  }
  if (getLevel(i_cells) > 0) {
    EG_BUG;
  }
  while (hasChildren(i_cells)) {
    vec3_t xc = getCellCentre(i_cells);
    if (x[0] > xc[0]) {
      if (x[1] > xc[1]) {
        if (x[2] > xc[2]) {
          i_cells = m_Cells[i_cells].m_Child[7];
        } else {
          i_cells = m_Cells[i_cells].m_Child[3];
        }
      } else {
        if (x[2] > xc[2]) {
          i_cells = m_Cells[i_cells].m_Child[5];
        } else {
          i_cells = m_Cells[i_cells].m_Child[1];
        }
      }
    } else {
      if (x[1] > xc[1]) {
        if (x[2] > xc[2]) {
          i_cells = m_Cells[i_cells].m_Child[6];
        } else {
          i_cells = m_Cells[i_cells].m_Child[2];
        }
      } else {
        if (x[2] > xc[2]) {
          i_cells = m_Cells[i_cells].m_Child[4];
        } else {
          i_cells = m_Cells[i_cells].m_Child[0];
        }
      }
    }
  }
  return i_cells;
}

bool Octree::intersectsFace(int cell, int face, vec3_t x1, vec3_t x2, double scale, double &k, double tol)
{
  vec3_t a, b, c;
  vec3_t x_centre = getCellCentre(cell);
  if (face == 0) {
    a = m_Nodes[m_Cells[cell].m_Node[0]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[2]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[4]].m_Position;
  } else if (face == 1) {
    a = m_Nodes[m_Cells[cell].m_Node[1]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[3]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[5]].m_Position;
  } else if (face == 2) {
    a = m_Nodes[m_Cells[cell].m_Node[0]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[1]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[4]].m_Position;
  } else if (face == 3) {
    a = m_Nodes[m_Cells[cell].m_Node[2]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[3]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[6]].m_Position;
  } else if (face == 4) {
    a = m_Nodes[m_Cells[cell].m_Node[0]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[1]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[2]].m_Position;
  } else if (face == 5) {
    a = m_Nodes[m_Cells[cell].m_Node[4]].m_Position;
    b = m_Nodes[m_Cells[cell].m_Node[5]].m_Position;
    c = m_Nodes[m_Cells[cell].m_Node[6]].m_Position;
  }
  a = x_centre + scale*(a - x_centre);
  b = x_centre + scale*(b - x_centre);
  c = x_centre + scale*(c - x_centre);
  vec3_t g1 = b-a;
  vec3_t g2 = c-a;
  double g1abs = g1.abs();
  double g2abs = g2.abs();
  g1.normalise();
  g2.normalise();
  vec3_t n = g1.cross(g2);
  k = GeometryTools::intersection(x1, x2-x1, a, n);
  bool intersects = false;
  if ((k > 0 - tol) && (k < 1 + tol)) {
    vec3_t x = x1 + k*(x2-x1) - a;
    double xg1 = x*g1;
    double xg2 = x*g2;
    xg1 /= g1abs;
    xg2 /= g2abs;
    if (fabs(x*n) > 1e-4) {
      EG_BUG;
    }
    if ((xg1 > 0 - tol) && (xg1 < 1 + tol) && (xg2 > 0 - tol) &&  (xg2 < 1 + tol)) {
      intersects = true;
    }
  }
  return intersects;
}

void Octree::resetRefineMarks()
{
  m_ToRefine.fill(false, m_Cells.size());
}

void Octree::getEdges(int cell, QVector<SortedPair<int> >& edges)
{
  edges.resize(12);
  edges[0].v1  = getNode(cell, 0); edges[0].v2  = getNode(cell, 1);
  edges[1].v1  = getNode(cell, 0); edges[1].v2  = getNode(cell, 2);
  edges[2].v1  = getNode(cell, 0); edges[2].v2  = getNode(cell, 4);
  edges[3].v1  = getNode(cell, 1); edges[3].v2  = getNode(cell, 3);
  edges[4].v1  = getNode(cell, 1); edges[4].v2  = getNode(cell, 5);
  edges[5].v1  = getNode(cell, 2); edges[5].v2  = getNode(cell, 3);
  edges[6].v1  = getNode(cell, 2); edges[6].v2  = getNode(cell, 6);
  edges[7].v1  = getNode(cell, 3); edges[7].v2  = getNode(cell, 7);
  edges[8].v1  = getNode(cell, 4); edges[8].v2  = getNode(cell, 5);
  edges[9].v1  = getNode(cell, 4); edges[9].v2  = getNode(cell, 6);
  edges[10].v1 = getNode(cell, 5); edges[10].v2 = getNode(cell, 7);
  edges[11].v1 = getNode(cell, 6); edges[11].v2 = getNode(cell, 7);
}

bool Octree::isInsideBounds(vec3_t x)
{
  bool inside = true;
  for (int i = 0; i < 3; ++i) {
    if ((x[i] < m_Corner1[i]) || (x[i] > m_Corner2[i])) {
      inside = false;
      break;
    }
  }
  return inside;
}

bool Octree::isInsideCell(int cell, vec3_t x, double overlap)
{
  vec3_t x0_1 = getNodePosition(cell, 0);
  vec3_t x0_2 = getNodePosition(cell, 7);
  vec3_t dx   = x0_2 - x0_1;
  vec3_t x1   = x0_1 - overlap*dx;
  vec3_t x2   = x0_2 + overlap*dx;

  for (int i = 0; i < 3; ++i) {
    if (x0_1[i] >= x0_2[i]) {
      EG_BUG;
    }
    if (x1[i] >= x2[i]) {
      EG_BUG;
    }
  }

  bool inside = true;
  for (int i = 0; i < 3; ++i) {
    if ((x[i] < x1[i]) || (x[i] > x2[i])) {
      inside = false;
      break;
    }
  }
  return inside;
}

bool Octree::triangleIntersectsCell(int cell, QVector<vec3_t> tri, double scale)
{
  if (tri.size() != 3) {
    EG_BUG;
  }

  // any node inside cell?
  foreach (vec3_t x, tri) {
    if (isInsideCell(cell, x, (scale - 1))) {
      return true;
    }
  }

  // any edge of cell intersects triangle?
  QVector<SortedPair<int> > edges;
  getEdges(cell, edges);
  foreach (SortedPair<int> edge, edges) {
    vec3_t x1 = getNodePosition(edge.v1);
    vec3_t x2 = getNodePosition(edge.v2);
    vec3_t xi, ri;
    vec3_t x_centre = getCellCentre(cell);
    x1 = x_centre + scale*(x1 - x_centre);
    x2 = x_centre + scale*(x2 - x_centre);
    if (GeometryTools::intersectEdgeAndTriangle(tri[0], tri[1], tri[2], x1, x2, xi, ri)) {
      return true;
    }
  }

  // any edge of triangle intersects any face of cell?
  for (int face = 0; face < 6; ++face) {
    double k;
    if (intersectsFace(cell, face, tri[0], tri[1], scale, k)) {
      return true;
    }
    if (intersectsFace(cell, face, tri[1], tri[2], scale, k)) {
      return true;
    }
    if (intersectsFace(cell, face, tri[2], tri[0], scale, k)) {
      return true;
    }
  }

  return false;
}

void Octree::getFinestChildren(int cell, QList<int> &finest_children)
{
  if (hasChildren(cell)) {
    for (int i = 0; i < 8; ++i) {
      getFinestChildren(m_Cells[cell].m_Child[i], finest_children);
    }
  } else {
    finest_children << cell;
  }
}

void Octree::getNeighbourRegion(int cell, QList<int> &neighbour_cells)
{
  QSet<int> cells;
  for (int i = 0; i < 8; ++i) {
    int node = getNode(cell, i);
    for (int j = 0; j < m_Node2Cell[node].size(); ++j) {
      int cell = m_Node2Cell[node][j];
      if (!hasChildren(cell)) {
        cells.insert(cell);
      }
    }
  }
  foreach (int cell, cells) {
    neighbour_cells.append(cell);
  }
}
