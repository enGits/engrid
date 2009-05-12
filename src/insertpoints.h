//
// C++ Interface: insertpoints
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
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

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class InsertPoints : public Operation
{
private:
  QMap <vtkIdType,bool> m_marked_cells;
  QMap <vtkIdType,bool> m_marked_nodes;
  QMap< pair<vtkIdType,vtkIdType>, vtkIdType> m_edge_map;
  
  QVector <stencil_t> m_StencilVector;
  QVector<vtkIdType> m_SelectedCells;
  QVector<vtkIdType> m_AllCells;
  QVector <vtkIdType> m_SelectedNodes;
  QVector <vtkIdType> m_AllNodes;
  
  QSet<int> m_bcs;
  
  int N_inserted_FP;
  int N_inserted_EP;
  int N_removed_FP;
  int N_removed_EP;
  
  int N_points;
  int N_cells;
  int N_newpoints;
  int N_newcells;
  int m_total_N_newpoints;
  int m_total_N_newcells;
  vtkIdType m_newNodeId;
  
  bool insert_FP;
  bool insert_EP;
  
/*  //for the UpdateDesiredMeshDensity operation
public:
  int MaxiterDensity;//used for UpdateDesiredMeshDensity operation
  void setMaxiterDensity(int a){MaxiterDensity=a;};
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void SetVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector){VMDvector=a_VMDvector;};*/
  
public:
  InsertPoints();

  ~InsertPoints();
  
  void operate();
  
  void Set_insert_FP(bool B){insert_FP=B;};
  void Set_insert_EP(bool B){insert_EP=B;};
  
  int insert_FP_counter();
  int insert_EP_counter();
  
  int insert_FP_actor(vtkUnstructuredGrid* grid_tmp);
  int insert_EP_actor(vtkUnstructuredGrid* grid_tmp);
  
  int insert_FP_all();
  int insert_EP_all();
  
  ///Check if a field point needs to be inserted
  bool insert_fieldpoint(vtkIdType D);
  ///Check if an edge point needs to be inserted
  bool insert_edgepoint(vtkIdType j,vtkIdType K);// node1 K, node2 j
};

#endif
