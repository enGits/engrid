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
#ifndef INSERTPOINTS_H
#define INSERTPOINTS_H

#include "surfaceoperation.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <QSet>
#include <QVector>
#include "egvtkobject.h"
#include "vertexmeshdensity.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

class InsertPoints : public SurfaceOperation
{
private:
  bool insert_FP;
  bool insert_EP;
  
  //attributes with setter functions
public:
  QSet<int> m_bcs;
  void setBCS(QSet<int> a_bcs) {m_bcs=a_bcs;};
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> const & a_VMDvector){VMDvector=a_VMDvector;};
  
public:
  InsertPoints();
  
  virtual void operate();
  
  void set_insert_FP(bool B){insert_FP=B;};
  void set_insert_EP(bool B){insert_EP=B;};
  
  int insert_FP_all();
  int insert_EP_all();
  
  ///Check if a field point needs to be inserted
  bool insert_fieldpoint(vtkIdType id_cell);
  
  ///Check if an edge point needs to be inserted
  bool insert_edgepoint(vtkIdType id_node1, vtkIdType id_node2);
  
  ///Check if an edge point needs to be inserted
  bool SplitSide(vtkIdType id_cell,int side);
  
  ///Returns the type of the node inserted on the edge S.p[1],S.p[3] from stencil_t S
  char getNewNodeType(stencil_t S);
  
};

#endif
