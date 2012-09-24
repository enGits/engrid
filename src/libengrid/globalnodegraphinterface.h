//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                      +
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

#ifndef GLOBALNODEGRAPHINTERFACE_H
#define GLOBALNODEGRAPHINTERFACE_H

#include "meshpartition.h"

/**
 * @brief An interface to run generic algorithms on the node graph of a grid.
 */
class GlobalNodeGraphInterface
{

  MeshPartition* m_Part;

public:

  typedef vtkIdType index_type;

  void   setMeshPartition(MeshPartition* part) { m_Part = part; }
  size_t size() { return m_Part->getGrid()->GetNumberOfPoints(); }
  size_t getNumLinks(size_t i) { return m_Part->n2nGSize(i); }
  size_t getLink(size_t i, size_t j) { return m_Part->n2nGG(i,j); }

};

#endif // GLOBALNODEGRAPHINTERFACE_H
