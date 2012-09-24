//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#ifndef SURFACENODEMOVEMENTCHECK_H
#define SURFACENODEMOVEMENTCHECK_H

#include "egvtkobject.h"
#include "meshpartition.h"

#include <QList>

class SurfaceNodeMovementCheck : public EgVtkObject
{

  vtkUnstructuredGrid* m_AuxGrid;
  vtkUnstructuredGrid* m_Grid;
  QVector<vtkIdType>   m_IdMap;
  MeshPartition        m_AuxPart;


public:

  SurfaceNodeMovementCheck();
  ~SurfaceNodeMovementCheck();

  void setGrid(vtkUnstructuredGrid* grid);
  bool checkNode(vtkIdType id_node, vec3_t x);
  void write(QString file_name);

};

#endif // SURFACENODEMOVEMENTCHECK_H
