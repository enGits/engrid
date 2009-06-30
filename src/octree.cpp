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
#include "octree.h"

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
}

void Octree::mergeNodes_identifyDuplicates()
{
  m_SameNodes.resize(m_Nodes.size());
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    m_SameNodes[i_nodes] = i_nodes;
  }
  foreach (const OctreeCell& cell, m_Cells) {
    double tol_cell = min(getDx(cell), min(getDy(cell), getDz(cell)));
    for (int i_neighbours = 0; i_neighbours < 6; ++i_neighbours) {
      if (cell.m_Neighbour[i_neighbours] != -1) {
        const OctreeCell& neigh = m_Cells[cell.m_Neighbour[i_neighbours]];
        double tol_neigh = min(getDx(neigh), min(getDy(neigh), getDz(neigh)));
        double tol = 0.01*min(tol_cell, tol_neigh);
        for (int i_nodes_cell = 0; i_nodes_cell < 8; ++i_nodes_cell) {
          for (int i_nodes_neigh = 0; i_nodes_neigh < 8; ++i_nodes_neigh) {
            int node_cell = cell.m_Node[i_nodes_cell];
            int node_neigh = neigh.m_Node[i_nodes_neigh];
            if (node_cell != node_neigh) {
              if ((m_Nodes[node_cell].m_Position - m_Nodes[node_neigh].m_Position).abs() < tol) {
                if (node_cell > node_neigh) {
                  m_SameNodes[node_cell] = node_neigh;
                } else {
                  m_SameNodes[node_neigh] = node_cell;
                }
              }
            }
          }
        }
      }
    }
  }

  QVector<bool> is_dup(m_Nodes.size(), false);
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    m_SameNodes[i_nodes] = m_SameNodes[m_SameNodes[i_nodes]];
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
      //cout << "DN: " << m_Nodes[i_nodes].m_Position << endl;
      ++N;
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, dup_grid);
  allocateGrid(dup_grid, N, N, false);
  N = 0;
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    if (is_dup[i_nodes]) {
      dup_grid->GetPoints()->SetPoint(N, m_Nodes[i_nodes].m_Position.data());
      vtkIdType pts[1];
      pts[0] = N;
      dup_grid->InsertNextCell(VTK_VERTEX, 1, pts);
      ++N;
    }
  }
  //writeGrid(dup_grid, "dup_nodes");
  cout << N << " duplicate nodes identified" << endl;
}

void Octree::mergeNodes_compactNodes()
{
  QVector<int> offset(m_Nodes.size());
  int last_offset = 0;
  for (int i = 0; i < m_Nodes.size(); ++i) {
    if (m_SameNodes[i] != i) {
      if (m_SameNodes[i] > i) {
        EG_BUG;
      }
      ++last_offset;
      //m_SameNodes[i] -= offset[m_SameNodes[i]];
    } else {
      //m_SameNodes[i] -= last_offset;
    }
    offset[i] = last_offset;
  }

  cout << "offset computed" << endl;

  for (int i = 0; i < m_Nodes.size(); ++i) {
    if (m_SameNodes[i] != i) {
      m_SameNodes[i] -= offset[m_SameNodes[i]];
    } else {
      m_SameNodes[i] -= offset[i];
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

void Octree::mergeNodes()
{
  //checkNeighbours();
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

void Octree::setBounds(vec3_t corner1, vec3_t corner2)
{
  if (m_Cells.size() != 1) {
    EG_BUG;
  }
  m_Corner1 = corner1;
  m_Corner2 = corner2;
  m_Dx = fabs(corner1[0] - corner2[0]);
  m_Dy = fabs(corner1[1] - corner2[1]);
  m_Dz = fabs(corner1[2] - corner2[2]);
  m_Nodes[0].m_Position = vec3_t(corner1[0], corner1[1], corner1[2]);
  m_Nodes[1].m_Position = vec3_t(corner2[0], corner1[1], corner1[2]);
  m_Nodes[2].m_Position = vec3_t(corner1[0], corner2[1], corner1[2]);
  m_Nodes[3].m_Position = vec3_t(corner2[0], corner2[1], corner1[2]);
  m_Nodes[4].m_Position = vec3_t(corner1[0], corner1[1], corner2[2]);
  m_Nodes[5].m_Position = vec3_t(corner2[0], corner1[1], corner2[2]);
  m_Nodes[6].m_Position = vec3_t(corner1[0], corner2[1], corner2[2]);
  m_Nodes[7].m_Position = vec3_t(corner2[0], corner2[1], corner2[2]);
  for (int i = 0; i < 8; ++i) {
    m_Cells[0].m_Node[i] = i;
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
  do {
    N1 = 0;
    N2 = 0;
    for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
      if (m_ToRefine[i_cells]) {
        ++N1;
        if (m_SmoothTransition) {
          for (int face = 0; face < 6; ++face) {
            //int neigh = getNeighbour(i_cells, face);
            int neigh = m_Cells[i_cells].m_Neighbour[face];
            if (neigh != -1) {
              if ((m_Cells[neigh].m_Level < m_Cells[i_cells].m_Level) && !m_ToRefine[neigh]) {
                m_ToRefine[neigh] = true;
                ++N2;
              }
            }
          }
        }
      }
    }
    ++count;
    if (count > 10) {
      cout << "ups!" << endl;
    }
  } while (N2 > 0);
  N2 = m_Cells.size();
  m_Cells.insert(N2, 8*N1, OctreeCell());
  int Nrefine = N1;
  int new_node = m_Nodes.size();
  m_Nodes.insert(m_Nodes.size(), N1*19, OctreeNode());
  N1 = N2;
  for (int i_cells = 0; i_cells < N1; ++i_cells) {
    if (m_ToRefine[i_cells]) {
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
  return Nrefine;
}

void Octree::toVtkGrid(vtkUnstructuredGrid *grid)
{
  int N = 0;
  for (int i = 0; i < m_Cells.size(); ++i) {
    if (m_Cells[i].m_Child[0] == -1) {
      ++N;
    }
  }
  allocateGrid(grid, N, m_Nodes.size());
  for (int id_node = 0; id_node < m_Nodes.size(); ++id_node) {
    grid->GetPoints()->SetPoint(id_node, m_Nodes[id_node].m_Position.data());
  }
  foreach (OctreeCell cell, m_Cells) {
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
        if (pts[i] < 0) {
          pts[i] = 0;
        }
      }
      grid->InsertNextCell(VTK_HEXAHEDRON, Npts, pts);
    }
  }
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
      EG_BUG;
    }
  }
  int i_cells = 0;
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

bool Octree::intersectsFace(int cell, int face, vec3_t x1, vec3_t x2, double tol)
{
  vec3_t a, b, c;
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
  vec3_t g1 = b-a;
  vec3_t g2 = c-a;
  double g1abs = g1.abs();
  double g2abs = g2.abs();
  g1.normalise();
  g2.normalise();
  vec3_t n = g1.cross(g2);
  double k = GeometryTools::intersection(x1, x2-x1, a, n);
  bool intersects = false;
  if ((k > 0 - tol*(x1-x2).abs()) && (k < 1 + tol*(x1-x2).abs())) {
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

