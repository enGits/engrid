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

  int    m_NumMaxIter;
  int    m_NumSmoothSteps;
  double m_MaxEdgeLength;
  double m_MinEdgeLength;
  double m_NodesPerQuarterCircle;
  double m_FeatureResolution2D;
  double m_FeatureResolution3D;
  bool   m_RespectFeatureEdgesForDeleteNodes;
  double m_FeatureAngleForDeleteNodes;
  bool   m_PerformGeometricTests;
  bool   m_UseProjectionForSmoothing;
  bool   m_UseNormalCorrectionForSmoothing;
  bool   m_AllowFeatureEdgeSwapping;
  double m_GrowthFactor;
  bool   m_SmoothSuccess;
  int    m_NumDelaunaySweeps;
  bool   m_AllowSmallAreaSwapping;
  bool   m_InsertNodes;
  bool   m_DeleteNodes;


protected: // methods

  void readSettings();


protected: // methods

  void prepare();
  void swap();
  void smooth(int N_iter, bool correct_curvature = false);
  int  insertNodes();
  int  deleteNodes();
  void computeMeshDensity();
  
public:

  SurfaceAlgorithm();

  void setVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector) { m_VMDvector = a_VMDvector; }
  void setMaxEdgeLength(double l)         { m_MaxEdgeLength = l; }
  void setNodesPerQuarterCircle(double N) { m_NodesPerQuarterCircle = N; }
  void setCellGrowthFactor(double cgf)    { m_GrowthFactor = cgf; }
  void setMaxNumIterations(int N)         { m_NumMaxIter = N; }
  void setNumSmoothSteps(int N)           { m_NumSmoothSteps = N; }
  void setNumDelaunaySweeps(int N)        { m_NumDelaunaySweeps = N; }
  void setDeleteNodesOn()                 { m_DeleteNodes = true; }
  void setDeleteNodesOff()                { m_DeleteNodes = false; }
  void setInsertNodesOn()                 { m_InsertNodes = true; }
  void setInsertNodesOff()                { m_InsertNodes = false; }
  void setDeleteNodes(bool s)             { m_DeleteNodes = s; }
  void setInsertNodes(bool s)             { m_InsertNodes = s; }

};

#endif // SURFACEALGORITHM_H
