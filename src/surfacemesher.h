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
#ifndef CREATESPECIALMAPPING_H
#define CREATESPECIALMAPPING_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCharArray.h>

#include <QSet>
#include <QVector>
#include <QString>
#include <QTextStream>
#include <QTime>

#include "egvtkobject.h"
#include "surfaceoperation.h"
#include "vertexmeshdensity.h"
#include "smoothingutilities.h"
#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"
#include "geometrytools.h"

#include <cmath>
#include <iostream>

class SurfaceMesher : public SurfaceOperation
{

private: // attributes

  QSet<int> m_BCs;
  QVector <VertexMeshDensity> VMDvector;
  int m_NumMaxIter;
  int m_NumSmoothSteps;

private: // methods

  void swap();
  void smooth(int N_iter);
  int insertNodes();
  int deleteNodes();

  void computeMeshDensity();
  void updateNodeInfo(bool update_type);

public:

  SurfaceMesher();

  virtual void operate();

  void setBoundaryCodes(QSet<int> bcs) { m_BCs = bcs; }

  //Used for UpdateDesiredMeshDensity operation
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector) { VMDvector = a_VMDvector; }

};

//end of SurfaceMesher class

#endif
