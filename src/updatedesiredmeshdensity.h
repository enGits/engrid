//
// C++ Interface: updatedesiredmeshdensity
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef UPDATEDESIREDMESHDENSITY_H
#define UPDATEDESIREDMESHDENSITY_H

#include "operation.h"

#include "vertexmeshdensity.h"

/// Update desired mesh density, i.e. the field used for surface meshing
class UpdateDesiredMeshDensity : public Operation
{
  //attributes
public:
  QSet<int> m_bcs;
  QVector <vtkIdType> m_SelectedNodes;
  QVector <vtkIdType> m_AllNodes;
  QVector<vtkIdType> m_SelectedCells;
  QVector<vtkIdType> m_AllCells;
  
  //attributes with setter functions
public:
  double Convergence_meshdensity;
  void SetConvergence_meshdensity(double C){Convergence_meshdensity=C;};
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void SetVertexMeshDensityVector(QVector <VertexMeshDensity> const & a_VMDvector){VMDvector=a_VMDvector;};
  int MaxiterDensity;//used for UpdateDesiredMeshDensity operation
  void setMaxiterDensity(int a){MaxiterDensity=a;};
  
  //methods
public:
    UpdateDesiredMeshDensity();

    ~UpdateDesiredMeshDensity();

  void operate();
  
  /// Get VertexMeshDensity object
  VertexMeshDensity getVMD(vtkIdType node, char VertexType);
};

#endif
