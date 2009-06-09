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

private: // attributes

  vec3_t m_Origin; /// origin of internal coordinate system
  mat3_t m_Base;   /// base vectors of internal coordinate system

public: // methods

  Octree(vec3_t origin = vec3_t(0,0,0),
         vec3_t basevec1 = vec3_t(1,0,0),
         vec3_t basevec2 = vec3_t(0,1,0),
         vec3_t basevec3 = vec3_t(0,0,1));

  vec3_t transfTo(vec3_t x);
  vec3_t transfFrom(vec3_t r);

};

#endif // OCTREE_H
