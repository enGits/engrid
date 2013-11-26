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
#ifndef UPDATEDESIREDMESHDENSITY_H
#define UPDATEDESIREDMESHDENSITY_H

#include "surfaceoperation.h"

#include "vertexmeshdensity.h"
#include "edgelengthsourcemanager.h"


/// Update desired mesh density, i.e. the field used for surface meshing

class UpdateDesiredMeshDensity : public SurfaceOperation
{

private: // data types

  struct point_t {
    vec3_t x;
    vec3_t n;
    double L;
    QList<int> idx;
    vtkIdType id_face;
  };

private: //attributes

  QSet<int>                   m_BCs;
  double                      m_GrowthFactor;
  double                      m_MaxEdgeLength;
  double                      m_MinEdgeLength;
  double                      m_NodesPerQuarterCircle;
  double                      m_FeatureResolution2D;
  double                      m_FeatureResolution3D;
  double                      m_FeatureThresholdAngle;
  QVector<double>             m_FeatureSize;
  QVector<bool>               m_Fixed;
  EdgeLengthSourceManager     m_ELSManager;
  bool                        m_OnlySurfaceCells;
  bool                        m_Relaxation;

protected: // methods

  void   computeFeature(const QList<point_t> points, QVector<double> &cl_pre, double res, int restriction_type);
  void   computeFeature2D(QVector<double> &cl_pre);
  void   computeFeature3D(QVector<double> &cl_pre);
  double computeSearchDistance(vtkIdType id_face);
  void   computeExistingLengths();


public: //methods

  UpdateDesiredMeshDensity();
  virtual void operate();

  void readSettings();
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> const & vmd) { m_VMDvector = vmd; }
  void setMaxEdgeLength(double l) { m_MaxEdgeLength = l; }
  void setMinEdgeLength(double l) { m_MinEdgeLength = l; }
  void setNodesPerQuarterCircle(double N) { m_NodesPerQuarterCircle = N; }
  void setCellGrowthFactor(double cgf) { m_GrowthFactor = cgf; }
  void setVolumeCellsOn() { m_OnlySurfaceCells = false; }
  void setVolumeCellsOff() { m_OnlySurfaceCells = true; }
  void setFeatureResolution2D(double n) { m_FeatureResolution2D = n; }
  void setFeatureResolution3D(double n) { m_FeatureResolution3D = n; }
  void setFeatureThresholdAngle(double a) { m_FeatureThresholdAngle = a; }
  void setRelaxationOff() { m_Relaxation = false; }
  void setRelaxationOn() { m_Relaxation = true; }

};

#endif
