// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "reducesurfacetriangulation.h"
#include "geometrytools.h"

ReduceSurfaceTriangulation::ReduceSurfaceTriangulation()
{
  EG_TYPENAME;
  m_PerformGeometricTests             = true;
  m_UseProjectionForSmoothing         = false;
  m_UseNormalCorrectionForSmoothing   = true;
  m_AllowFeatureEdgeSwapping          = true;
  m_RespectFeatureEdgesForDeleteNodes = true;
  m_FeatureAngleForDeleteNodes        = GeometryTools::deg2rad(20);

  m_NumDelaunaySweeps = 10;
  m_NumSmoothSteps = 1;
}

void ReduceSurfaceTriangulation::pass1()
{
  cout << "\nFirst pass of surface reduction:\n(This is the expensive part...)" << endl;
  int iter = 0;
  bool done = false;
  m_UseNormalCorrectionForSmoothing = true;
  int num_initial_nodes = m_Grid->GetNumberOfPoints();
  int num_del_max = 0;
  while (!done) {
    ++iter;
    cout << "\npass-1 iteration-" << iter << ":" << endl;
    computeMeshDensity();
    int num_deleted = deleteNodes();
    num_del_max = max(num_del_max, num_deleted);
    cout << "deleted nodes  : " << num_deleted << endl;
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
      smooth(1);
      swap();
    }
    done = num_deleted <= max(num_initial_nodes/100, num_del_max/100);
    cout << "total nodes : " << m_Grid->GetNumberOfPoints() << endl;
    cout << "total cells : " << m_Grid->GetNumberOfCells() << endl;
    //done = true;
  }
}

void ReduceSurfaceTriangulation::pass2()
{
  cout << "\n\nSecond pass of surface reduction:\n(This should be quick...)" << endl;
  m_UseNormalCorrectionForSmoothing = false;
  smooth(2);
}

void ReduceSurfaceTriangulation::operate()
{
  setStretchingFactor(1.0);
  prepare();
  //writeGrid(m_Grid, "take1");
  updateNodeInfo();
  //writeGrid(m_Grid, "take2");
  pass1();
  //pass2();
  createIndices(m_Grid);
  updateNodeInfo();
  computeMeshDensity();
}
