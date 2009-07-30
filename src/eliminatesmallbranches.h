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

#ifndef ELIMINATESMALLBRANCHES_H
#define ELIMINATESMALLBRANCHES_H

#include "operation.h"

class EliminateSmallBranches : public Operation
{

private: // attributes

  int             m_NumLayers;
  int             m_NumFillLayers;
  QVector<bool>   m_DeleteCell;
  QVector<bool>   m_IsSurfaceNode;
  QVector<bool>   m_MainVolumeCell;


protected: // methods

  bool needsToBeMarked(vtkIdType id_node, int layer = 0);
  void unmarkNode(vtkIdType id_node, int layer = 0);
  void fill(vtkIdType id_cell);
  void fillLayers();
  void fillCraters();

  virtual void operate();


public: // methods

  EliminateSmallBranches();

};

#endif // ELIMINATESMALLBRANCHES_H
