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
using namespace GeometryTools;

#include <cmath>
using namespace std;

#include <iostream>
using namespace std;

class SurfaceMesher : public SurfaceOperation {

private:
  int N_SmoothIterations;
  
  bool insert_FP;
  bool insert_EP;
  bool remove_FP;
  bool remove_EP;
  
  bool DoSwap;
  bool DoLaplaceSmoothing;
  
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
  QSet<int> m_bcs;
  
  double Convergence_meshdensity;
  
public:
  SurfaceMesher();
  virtual void operate();

  void setBoundaryCodes(QSet<int> a_bcs) { m_bcs=a_bcs; }

  //Used for UpdateDesiredMeshDensity operation
  int MaxiterDensity;//used for UpdateDesiredMeshDensity operation
  void setMaxiterDensity(int a) { MaxiterDensity=a; }
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  void setVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector) { VMDvector = a_VMDvector; }

  void setConvergence_meshdensity(double C) { Convergence_meshdensity = C; }
  void setDoSwap(bool B) { DoSwap = B; }
  void setDoLaplaceSmoothing (bool B) { DoLaplaceSmoothing = B; }
  void setN_SmoothIterations(int N) { N_SmoothIterations = N; }

  void set_insert_FP(bool B) { insert_FP = B; }
  void set_insert_EP(bool B) { insert_EP = B; }
  void set_remove_FP(bool B) { remove_FP = B; }
  void set_remove_EP(bool B) { remove_EP = B; }

  int SwapFunction();
  int SmoothFunction();
  void MeshDensityFunction();
  void UpdateNodeInfo(bool UpdateType);
};
//end of SurfaceMesher class

#endif
