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
#ifndef SURFACEALGORITHM_H
#define SURFACEALGORITHM_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCharArray.h>

#include <QSet>
#include <QVector>
#include <QString>
#include <QTextStream>
#include <QTime>

#include "surfaceoperation.h"
#include "vertexmeshdensity.h"
#include "geometrytools.h"

#include <cmath>
#include <iostream>

class SurfaceAlgorithm : public SurfaceOperation
{

protected: // attributes

  QVector <VertexMeshDensity> VMDvector;
  int    m_NumMaxIter;
  int    m_NumSmoothSteps;
  double m_MaxEdgeLength;
  double m_NodesPerQuarterCircle;
  bool   m_RespectFeatureEdgesForDeleteNodes;
  double m_FeatureAngleForDeleteNodes;
  bool   m_PerformGeometricTests;
  bool   m_UseProjectionForSmoothing;
  bool   m_UseNormalCorrectionForSmoothing;
  bool   m_AllowFeatureEdgeSwapping;
  double m_GrowthFactor;
  bool   m_SmoothSuccess;
  int    m_NumDelaunaySweeps;


private: // methods

  void readSettings();
  void readVMD();


protected: // methods

  void prepare();
  void swap();
  void smooth(int N_iter);
  int  insertNodes();
  int  deleteNodes();
  void computeMeshDensity();
  void updateNodeInfo(bool update_type = false);


public:

  SurfaceAlgorithm();

  void setVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector) { VMDvector = a_VMDvector; }
  void setMaxEdgeLength(double l)         { m_MaxEdgeLength = l; }
  void setNodesPerQuarterCircle(double N) { m_NodesPerQuarterCircle = N; }
  void setCellGrowthFactor(double cgf)    { m_GrowthFactor = cgf; }

};

#endif // SURFACEALGORITHM_H
