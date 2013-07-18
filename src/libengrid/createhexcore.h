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

#ifndef CREATEHEXCORE_H
#define CREATEHEXCORE_H

class CreateHexCore;

#include "operation.h"
#include "octree.h"

class CreateHexCore : public Operation
{

protected: // attributes

  vec3_t m_X1;
  vec3_t m_X2;
  vec3_t m_Xi;
  Octree m_Octree;
  int    m_NumI;
  int    m_NumJ;
  int    m_NumK;

protected: // methods

  virtual void operate();

  void refineOctree();
  void transferOctreeGrid();
  void deleteOutside(vtkUnstructuredGrid *grid);

public:

  CreateHexCore(vec3_t x1, vec3_t x2, vec3_t xi, int num_i, int num_j, int num_k);

};

#endif // CREATEHEXCORE_H
