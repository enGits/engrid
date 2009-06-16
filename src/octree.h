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
#ifndef OCTREE_H
#define OCTREE_H

#include "egvtkobject.h"

class Octree : public EgVtkObject
{

private: // data-types

  struct OctreeNode
  {
    vec3_t x;
  };

  struct OctreeCell
  {
    int n[8];
  };

private: // attributes

  vec3_t m_Origin;  /// origin of internal coordinate system
  mat3_t m_Base;    /// base vectors of internal coordinate system
  mat3_t m_InvBase; /// inverted base of internal coordiante system
  vec3_t m_Corner1; /// first corner of extend box of the whole domain (in internal coordinates)
  vec3_t m_Corner2; /// second corner of extend box of the whole domain (in internal coordinates)

  QVector<OctreeNode> m_Nodes;
  QVector<OctreeCell> m_Cells;

public: // methods

  Octree();

  void setOrigin(vec3_t x0);
  void setBase(vec3_t g1, vec3_t g2, vec3_t g3);
  void setBounds(vec3_t corner1, vec3_t corner2);

  vec3_t transfTo(vec3_t x);
  vec3_t transfFrom(vec3_t r);

};

inline vec3_t Octree::transfTo(vec3_t x)
{
  vec3_t dx = x - m_Origin;
  return m_InvBase*dx;
}

inline vec3_t Octree::transfFrom(vec3_t r)
{
  return m_Origin + m_Base*r;
}

#endif // OCTREE_H
