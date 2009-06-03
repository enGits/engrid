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
  void setConvergence_meshdensity(double C){Convergence_meshdensity=C;};
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> const & a_VMDvector){VMDvector=a_VMDvector;};
  int MaxiterDensity;//used for UpdateDesiredMeshDensity operation
  void setMaxiterDensity(int a){MaxiterDensity=a;};
  
  //methods
public:
    UpdateDesiredMeshDensity();

    ~UpdateDesiredMeshDensity();

    virtual void operate();
};

#endif
