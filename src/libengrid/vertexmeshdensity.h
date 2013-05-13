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
#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>
#include <QVector>
#include <QString>
#include "egvtkobject.h"

///\todo Rename density to edge length
class VertexMeshDensity
{

public: // methods

  VertexMeshDensity();/// Default constructor
  
  // node requirements
  /** Acceptable boundary codes for neighbour nodes. Unchecked=not allowed, checked=required
   * the first element is the boundary code
   * the second element is the acceptance value: 0=unchecked=not allowed, 1=partially checked=does not matter, 2=checked=required
   */
  QMap <int,int> BCmap;
  
  char type;/// Type of the node
  QSet <vtkIdType> nodeset;/// Set of acceptable node IDs
  // QVector<Qt::CheckState> BClist_value;// deprecated, replaced by BCmap
  
  // density
  double density;/// desired density for nodes matching all requirements
 
  // properties
  vtkIdType CurrentNode;/// ID of the current node
  
  /** operator to compare two VertexMeshDensity objects
   * this=user defined
   * VMD=VMD of current node
   * This operator is NOT SYMMETRICAL. But it can be used with indexOf.
 */
  bool operator==(const VertexMeshDensity & VMD) const;
  
  void setNodes(QString str);/// set nodeset by passing a string of the form "id1,id2,..."

  int findSmallestVMD( QVector <VertexMeshDensity> vector);
};

/// ostream operator to print out a VertexMeshDensity object
ostream& operator<<(ostream &out, VertexMeshDensity A);

/// ostream operator to print out a vector of VertexMeshDensity objects
ostream& operator<<(ostream &out, QVector<VertexMeshDensity> VMDvector);

#endif
