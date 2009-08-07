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

#include "surfacemesher.h"
#include "guimainwindow.h"

SurfaceMesher::SurfaceMesher() : SurfaceAlgorithm()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = true;
  m_UseNormalCorrectionForSmoothing = true;
  m_FeatureAngle = GeometryTools::deg2rad(200);
}

void SurfaceMesher::operate()
{
  prepare();
  EG_VTKDCN(vtkDoubleArray, md, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    md->SetValue(id_node, 1e-6);
  }
  updateNodeInfo(true);
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  while (!done) {
    ++iter;
    cout << "surface mesher iteration " << iter << ":" << endl;
    computeMeshDensity();
    num_inserted = insertNodes();
    cout << "  inserted nodes : " << num_inserted << endl;
    updateNodeInfo();
    swap();
    computeMeshDensity();
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
      smooth(1);
      swap();
    }
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    computeMeshDensity();
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
      smooth(1);
      swap();
    }
    //done = true;
    //done = (iter >= m_NumMaxIter);
    int N_crit = grid->GetNumberOfPoints()/100;
    done = (iter >= m_NumMaxIter) || ((num_inserted - num_deleted < N_crit) && (num_inserted + num_deleted < N_crit));
    //done = (iter >= m_NumMaxIter) || (num_inserted - num_deleted <= 0);
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  createIndices(grid);
  updateNodeInfo(false);
  computeMeshDensity();
  {
    int N1 = 0;
    int N2 = 0;
    QSet<int> bcs;
    GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
    foreach (int bc, bcs) {
      SurfaceProjection* proj = GuiMainWindow::pointer()->getSurfProj(bc);
      N1 += proj->getNumDirectProjections();
      N2 += proj->getNumFullSearches();
    }
    cout << N1 << " direct projections" << endl;
    cout << N2 << " full searches" << endl;
  }
}
