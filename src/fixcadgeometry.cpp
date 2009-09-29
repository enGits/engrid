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
//   setDesiredLength();
  mesher();

  // finalize
  createIndices(grid);
  updateNodeInfo(false);
  setDesiredLength();
}

void fixCadGeometry::mesher()
{
  setDesiredLength();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  while (!done) {
    ++iter;
    cout << "surface mesher iteration " << iter << ":" << endl;
    setDesiredLength();
    cout << "  inserted nodes : " << num_inserted << endl;
    updateNodeInfo();
    setDesiredLength();
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    setDesiredLength();
    int N_crit = grid->GetNumberOfPoints()/100;
    done = (iter >= m_NumMaxIter) || ((num_inserted - num_deleted < N_crit) && (num_inserted + num_deleted < N_crit));
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  createIndices(grid);
  updateNodeInfo(false);
  setDesiredLength();
}

void fixCadGeometry::setDesiredLength(double L)
{
  setAllSurfaceCells();
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2n   = getPartN2N();
  l2l_t  c2c   = getPartC2C();
  
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    characteristic_length_specified, grid, "node_specified_density");
  
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    characteristic_length_specified->SetValue(id_node, 0);
    characteristic_length_desired->SetValue(id_node, L);
  }
}
