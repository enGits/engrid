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
#ifndef CREATEHEXIBMESH_H
#define CREATEHEXIBMESH_H

class CreateHexIbMesh;

#include "surfaceoperation.h"
#include "octree.h"
#include "edgelengthsourcemanager.h"

class CreateHexIbMesh : public SurfaceOperation
{

  int                        m_MinDim;
  double                     m_MinSize;
  Octree                     m_Octree;
  QVector<QList<vtkIdType> > m_Faces;
  QVector<double>            m_MeshSize;
  int                        m_MinNumLayersWithRequiredResolution;
  vec3_t                     m_InsidePosition;
  EdgeLengthSourceManager    m_ELSManager;
  double                     m_GrowthFactor;
  double                     m_MinEdgeLength;
  double                     m_MaxEdgeLength;


protected: // methods

  int    refine();
  void   updateMeshSize();
  double meshSize(vtkIdType id_face);
  double meshSize(const QList<vtkIdType> &faces);
  void   findInsideCells(MeshPartition &part, QList<vtkIdType> &inside_cells);

  virtual void operate();

  QString bigIntText(long long N);


public:

  CreateHexIbMesh();
  void setMinNumLayersWithRequiredResolution(int N) { m_MinNumLayersWithRequiredResolution = N; }
  void setMinDim(int N) { m_MinDim = N; }
  void setInsidePosition(vec3_t x) { m_InsidePosition = x; }

};

#endif // CREATEHEXIBMESH_H
