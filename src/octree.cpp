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

Octree::Octree()
{
  m_Origin = vec3_t(0,0,0);
  setBase(vec3_t(1,0,0), vec3_t(0,1,0), vec3_t(0,0,1));
  m_Nodes.resize(8);
  m_Cells.resize(1);
  setBounds(vec3_t(0,0,0), vec3_t(1,1,1));
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
}


