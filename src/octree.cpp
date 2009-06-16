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




Octree::Octree()
{
  m_Origin = vec3_t(0,0,0);
  setBase(vec3_t(1,0,0), vec3_t(0,1,0), vec3_t(0,0,1));
  m_Nodes.resize(8);
  m_Cells.resize(1);
  setBounds(vec3_t(0,0,0), vec3_t(1,1,1));
  setSmoothTransitionOn();
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
    for (int cell = 0; cell < m_Cells.size(); ++cell) {
      if (m_ToRefine[cell]) {
        ++N1;
        for (int face = 0; face < 6; ++face) {
          int neigh = getNeighbour(cell, face);
          if (neigh != -1) {
            if (m_Cells[neigh].m_Level < m_Cells[cell].m_Level) {
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
  N1 = N2;
  for (int cell = 0; cell < N1; ++cell) {
    if (m_ToRefine[cell]) {
      for (int child = 0; child < 8; ++child) {
        m_Cells[cell].m_Child[child] = N2;
        m_Cells[N2].m_Parent = cell;
        m_Cells[N2].m_Level = m_Cells[cell].m_Level + 1;
        ++N2;
      }
      // - - -
      m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[1] = m_Cells[cell].m_Child[1];
      m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[3] = m_Cells[cell].m_Child[2];
      m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[5] = m_Cells[cell].m_Child[4];
      // + - -
      m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[0] = m_Cells[cell].m_Child[0];
      m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[3] = m_Cells[cell].m_Child[3];
      m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[5] = m_Cells[cell].m_Child[5];
      // - + -
      m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[1] = m_Cells[cell].m_Child[3];
      m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[2] = m_Cells[cell].m_Child[0];
      m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[5] = m_Cells[cell].m_Child[6];
      // + + -
      m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[0] = m_Cells[cell].m_Child[2];
      m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[2] = m_Cells[cell].m_Child[1];
      m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[5] = m_Cells[cell].m_Child[7];
      // - - +
      m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[1] = m_Cells[cell].m_Child[5];
      m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[3] = m_Cells[cell].m_Child[6];
      m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[4] = m_Cells[cell].m_Child[0];
      // + - +
      m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[0] = m_Cells[cell].m_Child[4];
      m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[3] = m_Cells[cell].m_Child[7];
      m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[4] = m_Cells[cell].m_Child[1];
      // - + +
      m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[1] = m_Cells[cell].m_Child[7];
      m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[2] = m_Cells[cell].m_Child[4];
      m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[4] = m_Cells[cell].m_Child[2];
      // + + +
      m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[0] = m_Cells[cell].m_Child[6];
      m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[2] = m_Cells[cell].m_Child[5];
      m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[4] = m_Cells[cell].m_Child[3];
    }
  }
  for (int cell = 0; cell < N1; ++cell) {
    if (m_Cells[cell].m_Child[0] >= N1) {
      {
        int neigh = m_Cells[cell].m_Neighbour[0];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[0] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[0] = m_Cells[neigh].m_Child[3];
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[0] = m_Cells[neigh].m_Child[5];
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[0] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[0] = neigh;
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[0] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[cell].m_Neighbour[1];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[1] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[1] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[1] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[1] = m_Cells[neigh].m_Child[6];
          } else {
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[1] = neigh;
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[1] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[cell].m_Neighbour[2];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[2] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[2] = m_Cells[neigh].m_Child[3];
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[2] = m_Cells[neigh].m_Child[6];
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[2] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[2] = neigh;
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[2] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[cell].m_Neighbour[3];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[3] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[3] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[3] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[3] = m_Cells[neigh].m_Child[5];
          } else {
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[3] = neigh;
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[3] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[cell].m_Neighbour[4];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[4] = m_Cells[neigh].m_Child[4];
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[4] = m_Cells[neigh].m_Child[5];
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[4] = m_Cells[neigh].m_Child[6];
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[4] = m_Cells[neigh].m_Child[7];
          } else {
            m_Cells[m_Cells[cell].m_Child[0]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[cell].m_Child[1]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[cell].m_Child[2]].m_Neighbour[4] = neigh;
            m_Cells[m_Cells[cell].m_Child[3]].m_Neighbour[4] = neigh;
          }
        }
      }
      {
        int neigh = m_Cells[cell].m_Neighbour[5];
        if (neigh != -1) {
          if (m_Cells[neigh].m_Child[0] != -1) {
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[5] = m_Cells[neigh].m_Child[0];
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[5] = m_Cells[neigh].m_Child[1];
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[5] = m_Cells[neigh].m_Child[2];
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[5] = m_Cells[neigh].m_Child[3];
          } else {
            m_Cells[m_Cells[cell].m_Child[4]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[cell].m_Child[5]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[cell].m_Child[6]].m_Neighbour[5] = neigh;
            m_Cells[m_Cells[cell].m_Child[7]].m_Neighbour[5] = neigh;
          }
        }
      }
    }
  }
}
