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
  setBounds(vec3_t(0,0,0), vec3_t(1,1,1));
  setSmoothTransitionOn();
}

void Octree::resetNodeMerge()
{
  m_SameNodes.resize(m_Nodes.size());
  for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
    m_SameNodes[i_nodes] = i_nodes;
  }
}

void Octree::mergeNodes()
{
  foreach (const OctreeCell& cell, m_Cells) {
    for (int i_neighbours = 0; i_neighbours < 6; ++i_neighbours) {
      const OctreeCell& neigh = m_Cells[cell.m_Neighbour[i_neighbours]];
      for (int i_nodes_cell = 0; i_nodes_cell < 8; ++i_nodes_cell) {
        for (int i_nodes_neigh = 0; i_nodes_neigh < 8; ++i_nodes_neigh) {
        }
      }
    }
  }
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

void Octree::refineAll()
{
  int N1;
  int N2;
  do {
    N1 = 0;
    N2 = 0;
    for (int i_cells = 0; i_cells < m_Cells.size(); ++i_cells) {
      if (m_ToRefine[i_cells]) {
        ++N1;
        for (int face = 0; face < 6; ++face) {
          int neigh = getNeighbour(i_cells, face);
          if (neigh != -1) {
            if (m_Cells[neigh].m_Level < m_Cells[i_cells].m_Level) {
              m_ToRefine[neigh] = true;
              ++N2;
            }
          }
        }
      }
    }
  } while (N2 > 0);
  N2 = m_Cells.size();
  m_Cells.insert(N2, N1, OctreeCell());
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
      ne[0] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position);
      ne[1] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position);
      ne[2] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position);
      ne[3] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position);
      ne[4] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position);
      ne[5] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      ne[6] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position);
      ne[7] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      ne[8] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[0]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[4]].m_Position);
      ne[9] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[1]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[5]].m_Position);
      ne[10] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[2]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[6]].m_Position);
      ne[11] = new_node++;
      m_Nodes[new_node].m_Position = 0.5*(m_Nodes[m_Cells[i_cells].m_Node[3]].m_Position + m_Nodes[m_Cells[i_cells].m_Node[7]].m_Position);
      int nf[6];
      nf[0] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[0]].m_Position + m_Nodes[ne[4]].m_Position
                                           + m_Nodes[ne[8]].m_Position  + m_Nodes[ne[10]].m_Position);
      nf[1] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[1]].m_Position + m_Nodes[ne[5]].m_Position
                                           + m_Nodes[ne[9]].m_Position  + m_Nodes[ne[11]].m_Position);
      nf[2] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[2]].m_Position + m_Nodes[ne[6]].m_Position
                                           + m_Nodes[ne[8]].m_Position  + m_Nodes[ne[9]].m_Position);
      nf[3] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[3]].m_Position + m_Nodes[ne[7]].m_Position
                                           + m_Nodes[ne[10]].m_Position + m_Nodes[ne[11]].m_Position);
      nf[4] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[0]].m_Position + m_Nodes[ne[1]].m_Position
                                           + m_Nodes[ne[2]].m_Position  + m_Nodes[ne[3]].m_Position);
      nf[5] = new_node++;
      m_Nodes[new_node].m_Position = 0.25*(  m_Nodes[ne[4]].m_Position + m_Nodes[ne[5]].m_Position
                                           + m_Nodes[ne[6]].m_Position  + m_Nodes[ne[7]].m_Position);
      int nv = new_node++;
      m_Nodes[new_node].m_Position = 1.0/6.0*(  m_Nodes[nf[0]].m_Position + m_Nodes[nf[1]].m_Position + m_Nodes[nf[2]].m_Position
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
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[0] = nf[0];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[1] = nv;
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[2] = ne[10];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[3] = nf[3];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[4] = ne[4];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[5] = nf[5];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[6] = nn[6];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[7] = ne[7];
      // child 7
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[0] = nv;
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[1] = nf[1];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[2] = nf[3];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[3] = ne[11];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[4] = nf[5];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[5] = ne[5];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[6] = ne[7];
      m_Cells[m_Cells[i_cells].m_Child[0]].m_Node[7] = nn[7];

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
          if (m_Cells[neigh].m_Child[0] != -1) {
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
}
