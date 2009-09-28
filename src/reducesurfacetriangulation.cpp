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

#include "reducesurfacetriangulation.h"

ReduceSurfaceTriangulation::ReduceSurfaceTriangulation()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = false;
  m_UseNormalCorrectionForSmoothing = true;
  m_FeatureAngle = GeometryTools::deg2rad(15);
  m_AllowFeatureEdgeSwapping = false;
}

void ReduceSurfaceTriangulation::pass1()
{
  cout << "\nFirst pass of surface reduction:\n(This is the expensive part...)" << endl;
  int iter = 0;
  bool done = false;
  m_UseNormalCorrectionForSmoothing = true;
  m_RespectFeatureEdgesForDeleteNodes = false;
  int num_initial_nodes = grid->GetNumberOfPoints();
  int num_del_max = 0;
  while (!done) {
    ++iter;
    cout << "\npass-1 iteration-" << iter << ":" << endl;
    cout << "computing 'snap-points'" << endl;
    UpdatePotentialSnapPoints(true, false);
    cout << "computing characteristic length" << endl;
    computeMeshDensity();
    cout << "removing nodes" << endl;
    int num_deleted = deleteNodes();
    num_del_max = max(num_del_max, num_deleted);
    cout << "deleted nodes  : " << num_deleted << endl;
    cout << "performing delaunay swap" << endl;
    swap();
    cout << "computing 'snap-points'" << endl;
    UpdatePotentialSnapPoints(true, false);
    cout << "smoothing" << endl;
    smooth(1);
    done = num_deleted <= num_del_max/100;
    cout << "total nodes : " << grid->GetNumberOfPoints() << endl;
    cout << "total cells : " << grid->GetNumberOfCells() << endl;
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
  prepare();
  writeGrid(grid, "take1");
  updateNodeInfo(true);
  writeGrid(grid, "take2");
  pass1();
  pass2();
  createIndices(grid);
  updateNodeInfo(false);
  computeMeshDensity();
}
