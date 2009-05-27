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
#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>
#include <QVector>
#include <QString>
#include "egvtkobject.h"

class VertexMeshDensity{
public:
  VertexMeshDensity();
//   ~VertexMeshDensity();

public:
  QVector <Qt::CheckState> BClist_value;
  QMap <int,int> BCmap;
  char type;
  QSet <vtkIdType> nodeset;
  double density;
  vtkIdType CurrentNode;
  bool operator==(const VertexMeshDensity & VMD) const;
  void SetNodes(QString str);
};

ostream& operator<<(ostream &out, VertexMeshDensity A);

ostream& operator<<(ostream &out, QVector<VertexMeshDensity> VMDvector);

#endif
