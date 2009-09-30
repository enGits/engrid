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

#include <vtkMath.h>

#include "geometrytools.h"
using namespace GeometryTools;

#include <QtDebug>
#include <iostream>
using namespace std;

FixCadGeometry::FixCadGeometry()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = false;
  m_UseNormalCorrectionForSmoothing = true;
  m_FeatureAngle = GeometryTools::deg2rad(0.5);
  m_AllowFeatureEdgeSwapping = false;
  m_AllowSmallAreaSwapping = true;
  m_NumDelaunaySweeps = 100;
}

void FixCadGeometry::operate()
{
  setAllCells();

  //prepare BCmap
  QSet <int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  QMap <int,int> BCmap;
  foreach(int bc, bcs) BCmap[bc]=1;

  //set density infinite
  VertexMeshDensity VMD;
  VMD.density = 1e99;
  VMD.BCmap = BCmap;
  cout<<"VMD="<<VMD<<endl;
  qWarning()<<"VMD.BCmap="<<VMD.BCmap;
  m_VMDvector.push_back(VMD);
  cout<<"m_VMDvector="<<m_VMDvector<<endl;

  customUpdateNodeInfo();

  //call surface mesher
  setGrid(grid);
  setBoundaryCodes(bcs);
  setVertexMeshDensityVector(m_VMDvector);
  setDesiredLength();
  mesher();

  // finalise
  createIndices(grid);
  customUpdateNodeInfo();
  setDesiredLength();
}

void FixCadGeometry::mesher()
{
  setDesiredLength();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  swap();
  while (!done) {
    ++iter;
    cout << "CAD fix iteration " << iter << ":" << endl;
    setDesiredLength();
    customUpdateNodeInfo();
    setDesiredLength();
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    setDesiredLength();
    done = (iter >= 20) || (num_deleted == 0);
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  createIndices(grid);
  customUpdateNodeInfo();
  setDesiredLength();
}

void FixCadGeometry::setDesiredLength(double L)
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

void FixCadGeometry::customUpdateNodeInfo()
{
  setAllCells();
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  l2g_t cells = getPartCells();
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    node_type->SetValue(id_node, VTK_FIXED_VERTEX);
  }
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, grid)) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      if (N_pts == 3) {
        vec3_t n1 = cellNormal(grid, id_cell);
        QVector<int> num_bad_edges(3, 0);
        for (int i = 0; i < N_pts; ++i) {
          int i1 = i;
          int i2 = 0;
          if (i < N_pts - 1) {
            i2= i+1;
          }
          vec3_t x1, x2;
          grid->GetPoint(pts[i1], x1.data());
          grid->GetPoint(pts[i2], x2.data());
          double L = (x1 - x2).abs();
          vtkIdType id_neigh_cell = m_Part.c2cGG(id_cell, i);
          if (id_neigh_cell != -1) {
            bool bad_edge = false;
            vec3_t n2 = cellNormal(grid, id_neigh_cell);
            if (angle(n1, n2) > deg2rad(179.5)) {
              bad_edge = true;
            }
            if (bad_edge) {
              ++num_bad_edges[i1];
              ++num_bad_edges[i2];
            }
          }
        }
        for (int i = 0; i < N_pts; ++i) {
          if (num_bad_edges[i] >= 2) {
            node_type->SetValue(pts[i], VTK_SIMPLE_VERTEX);
          }
        }
      }
    }
  }
}

