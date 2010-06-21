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

#include "surfacemesher.h"
#include "guimainwindow.h"

#include "laplacesmoother.h"

SurfaceMesher::SurfaceMesher() : SurfaceAlgorithm()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  //m_UseProjectionForSmoothing = true;
  getSet("surface meshing", "use surface projection for smoothing", true, m_UseProjectionForSmoothing);
  getSet("surface meshing", "use normal correction for smoothing", false, m_UseNormalCorrectionForSmoothing);
  getSet("surface meshing", "allow feature edge swapping", false, m_AllowFeatureEdgeSwapping);
  //m_UseNormalCorrectionForSmoothing = true;
  //m_AllowFeatureEdgeSwapping = false;
  m_EdgeAngle = m_FeatureAngle;
  
  getSet("surface meshing", "interpolate after meshing (experimental)", false, m_interpolateAfterMeshing);
}

void SurfaceMesher::operate()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  computeMeshDensity();
  prepare();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
  updateNodeInfo(true);
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = (iter >= m_NumMaxIter);
  int Nfull = 0;
  int Nhalf = 0;
  while (!done) {
    ++iter;
    cout << "surface mesher iteration " << iter << ":" << endl;
    computeMeshDensity();
    //return;
    num_inserted = insertNodes();
    cout << "  inserted nodes : " << num_inserted << endl;
    updateNodeInfo();
    swap();
    computeMeshDensity();
    num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    computeMeshDensity();
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
      SurfaceProjection::Nfull = 0;
      SurfaceProjection::Nhalf = 0;
      smooth(1);
      Nfull += SurfaceProjection::Nfull;
      Nhalf += SurfaceProjection::Nhalf;
      cout << "    " << SurfaceProjection::Nfull << " full searches" << endl;
      cout << "    " << SurfaceProjection::Nhalf << " half searches" << endl;
      swap();
    }
    int N_crit = m_Grid->GetNumberOfPoints()/100;
    //done = (iter >= m_NumMaxIter) || ((num_inserted - num_deleted < N_crit) && (num_inserted + num_deleted < N_crit));
    done = (iter >= m_NumMaxIter);
    cout << "  total nodes    : " << m_Grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << m_Grid->GetNumberOfCells() << endl;
  }
  createIndices(m_Grid);
  updateNodeInfo(false);
  computeMeshDensity();
  {
    QSet<int> bcs;
    GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
    foreach (int bc, bcs) {
      SurfaceProjection* proj = GuiMainWindow::pointer()->getSurfProj(bc);
    }
  }
  cout << Nfull << " full searches in total" << endl;
  cout << Nhalf << " half searches in total" << endl;

  if(m_interpolateAfterMeshing) {
    qDebug()<<"+++ CORRECTING CURVATURE +++";
    // correct curvature
    LaplaceSmoother lap;
    lap.setCorrectCurvature(true);
    lap.setNoCheck(true);
    lap.setGrid(m_Grid);
    QVector<vtkIdType> cls;
    getSurfaceCells(m_BoundaryCodes, cls, m_Grid);
    lap.setCells(cls);
    lap.setNumberOfIterations(1);
    lap.setProjectionIterations(2);
    m_BoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();
    lap.setBoundaryCodes(m_BoundaryCodes);//IMPORTANT: so that unselected nodes become fixed when node types are updated!
    if (m_UseProjectionForSmoothing) {
      lap.setProjectionOn();
    } else {
      lap.setProjectionOff();
    }
    if (m_UseNormalCorrectionForSmoothing) {
      lap.setNormalCorrectionOn();
    } else {
      lap.setNormalCorrectionOff();
    }
    lap();
    m_SmoothSuccess = lap.succeeded();
  }
}
