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
  m_RespectFeatureEdgesForDeleteNodes = true;
  m_FeatureAngleForDeleteNodes = GeometryTools::deg2rad(90.0);
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = false;
}

void ReduceSurfaceTriangulation::operate()
{
  prepare();
  updateNodeInfo(true);
  int iter = 0;
  bool done = false;
  while (!done) {
    ++iter;
    cout << "reduce surface triangulation iteration " << iter << ":" << endl;
    computeMeshDensity();
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    //smooth(1);
    computeMeshDensity();
    done = num_deleted == 0;
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  /*
  for (int i = 0; i < 3; ++i) {
    smooth(1);
    swap();
  }
  */
  createIndices(grid);
  updateNodeInfo(false);
  computeMeshDensity();
  writeGrid(grid, "after_reduction");
}
