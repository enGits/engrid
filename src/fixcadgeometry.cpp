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
#include "fixcadgeometry.h"
#include "surfacemesher.h"
#include "vertexmeshdensity.h"
#include "guimainwindow.h"

#include <QtDebug>
#include <iostream>
using namespace std;

fixCadGeometry::fixCadGeometry()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = true;
  m_UseNormalCorrectionForSmoothing = true;
  m_FeatureAngle = GeometryTools::deg2rad(200);//this angle is also used by swaptriangles!!!
}

void fixCadGeometry::operate()
{
  qDebug()<<"==>fixing CAD geometry...";
//   prepare();
  setAllCells();
  readSettings();
//   readVMD();
  
  //prepare BCmap
  QSet <int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  QMap <int,int> BCmap;
  foreach(int bc, bcs) BCmap[bc]=1;
  
  //set density infinite
  VertexMeshDensity VMD;
  VMD.density = 9000;
  VMD.BCmap = BCmap;
  cout<<"VMD="<<VMD<<endl;
  qWarning()<<"VMD.BCmap="<<VMD.BCmap;
  m_VMDvector.push_back(VMD);
  cout<<"m_VMDvector="<<m_VMDvector<<endl;
  
  //update node info
  updateNodeInfo(true);
  
  //call surface mesher
  setGrid(grid);
  setBoundaryCodes(bcs);
  setVertexMeshDensityVector(m_VMDvector);
  computeMeshDensity();
//   mesher();

  // finalize
  createIndices(grid);
  updateNodeInfo(false);
  computeMeshDensity();
}

void fixCadGeometry::mesher()
{
/*  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }*/
  computeMeshDensity();
//   prepare();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
//   updateNodeInfo(true);
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  while (!done) {
    ++iter;
    cout << "surface mesher iteration " << iter << ":" << endl;
    computeMeshDensity();
//     num_inserted = insertNodes();
    cout << "  inserted nodes : " << num_inserted << endl;
    updateNodeInfo();
//     swap();
    computeMeshDensity();
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
//       smooth(1);
//       swap();
    }
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    computeMeshDensity();
    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      cout << "  smoothing    : " << i+1 << "/" << m_NumSmoothSteps << endl;
//       smooth(1);
//       swap();
    }
    int N_crit = grid->GetNumberOfPoints()/100;
    done = (iter >= m_NumMaxIter) || ((num_inserted - num_deleted < N_crit) && (num_inserted + num_deleted < N_crit));
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
