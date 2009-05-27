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

#include "operation.h"

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

class InsertPoints : public Operation
{
private:
  QMap <vtkIdType,bool> m_marked_cells;
  
  QVector <stencil_t> m_StencilVector;
  QVector<vtkIdType> m_SelectedCells;
  QVector<vtkIdType> m_AllCells;
  QVector <vtkIdType> m_SelectedNodes;
  QVector <vtkIdType> m_AllNodes;
  
  int m_N_points;
  int m_N_cells;
  int m_N_newpoints;
  int m_N_newcells;
  
  int m_total_N_newpoints;
  int m_total_N_newcells;
  vtkIdType m_newNodeId;
  
  bool insert_FP;
  bool insert_EP;
  
  //attributes with setter functions
public:
  QSet<int> m_bcs;
  void SetBCS(QSet<int> a_bcs) {m_bcs=a_bcs;};
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void SetVertexMeshDensityVector(QVector <VertexMeshDensity> const & a_VMDvector){VMDvector=a_VMDvector;};
  
public:
  InsertPoints();
  
  void operate();
  
  void Set_insert_FP(bool B){insert_FP=B;};
  void Set_insert_EP(bool B){insert_EP=B;};
  
  int insert_FP_counter();
  int insert_FP_actor(vtkUnstructuredGrid* grid_tmp);
  int insert_FP_all();
  
  int insert_EP_all();
  
  ///Check if a field point needs to be inserted
  bool insert_fieldpoint(vtkIdType D);
  ///Check if an edge point needs to be inserted
  bool insert_edgepoint(vtkIdType j,vtkIdType K);// node1 K, node2 j

  double NewCurrentMeshDensity(vtkIdType a_vertex,double a_dist);

};

#endif
