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

class Octree;

#include "egvtkobject.h"

class OctreeNode
{

  friend class Octree;

  vec3_t m_Position;

public:

  vec3_t& getPosition()                { return m_Position; }
  void    setPosition(const vec3_t& x) { m_Position = x; }

};

class OctreeCell
{

  friend class Octree;

  int m_Node[8];
  int m_Child[8];
  int m_Neighbour[6];
  int m_Parent;
  int m_Level;

public:

  OctreeCell();

  int  getNode     (int i) { return m_Node[i]; }
  int  getNeighbour(int i) { return m_Neighbour[i]; }

};

class Octree : public EgVtkObject
{

private: // attributes

  vec3_t m_Origin;  ///< origin of internal coordinate system
  mat3_t m_Base;    ///< base vectors of internal coordinate system
  mat3_t m_InvBase; ///< inverted base of internal coordiante system
  vec3_t m_Corner1; ///< first corner of extend box of the whole domain (in internal coordinates)
  vec3_t m_Corner2; ///< second corner of extend box of the whole domain (in internal coordinates)
  bool   m_SmoothTransition;

  QVector<OctreeNode> m_Nodes;
  QVector<OctreeCell> m_Cells;
  QVector<bool>       m_ToRefine;

public: // methods

  Octree();

  void setOrigin(vec3_t x0);
  void setBase(vec3_t g1, vec3_t g2, vec3_t g3);
  void setBounds(vec3_t corner1, vec3_t corner2);

  vec3_t transfTo(vec3_t x);
  vec3_t transfFrom(vec3_t r);

  int  getNeighbour(int cell, int neigh) { return m_Cells[cell].m_Neighbour[neigh]; }

  void markToRefine(int cell) { m_ToRefine[cell] = true; }
  void refineAll();
  void setSmoothTransitionOn()  { m_SmoothTransition = true; }
  void setSmoothTransitionOff() { m_SmoothTransition = false; }

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
