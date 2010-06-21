//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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

private: //attributes

  QSet<int>                   m_BCs;
  double                      m_GrowthFactor;
  QVector <VertexMeshDensity> m_VMDvector; ///< the mesh density rules
  double                      m_MaxEdgeLength;
  double                      m_NodesPerQuarterCircle;
  int                         m_MinMumCellsAcross;
  QVector<bool>               m_Fixed;
  EdgeLengthSourceManager     m_ELSManager;


protected: // methods

  void computeExistingLengths();


public: //methods

  UpdateDesiredMeshDensity();
  virtual void operate();
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> const & vmd) { m_VMDvector = vmd; }
  void setMaxEdgeLength(double l) { m_MaxEdgeLength = l; }
  void setNodesPerQuarterCircle(double N) { m_NodesPerQuarterCircle = N; }
  void setCellGrowthFactor(double cgf) { m_GrowthFactor = cgf; }

};

#endif
