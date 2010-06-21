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
#include "fixcadgeometry.h"
#include "surfacemesher.h"
#include "vertexmeshdensity.h"
#include "guimainwindow.h"
#include "deletestraynodes.h"
#include "correctsurfaceorientation.h"

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

void FixCadGeometry::callMesher()
{
  setDesiredLength();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
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
    cout << "  total nodes    : " << m_Grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << m_Grid->GetNumberOfCells() << endl;
  }
  createIndices(m_Grid);
  customUpdateNodeInfo();
  setDesiredLength();
}

void FixCadGeometry::setDesiredLength(double L)
{
  setAllSurfaceCells();
  l2g_t  nodes = getPartNodes();

  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   m_Grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    characteristic_length_specified, m_Grid, "node_specified_density");

  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    characteristic_length_specified->SetValue(id_node, 0);
    characteristic_length_desired->SetValue(id_node, L);
  }
}

void FixCadGeometry::customUpdateNodeInfo()
{
  cout << "updating node information ..." << endl;
  setAllCells();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  l2g_t cells = getPartCells();
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    node_type->SetValue(id_node, VTK_FIXED_VERTEX);
  }
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      if (N_pts == 3) {
        vec3_t n1 = cellNormal(m_Grid, id_cell);
        QVector<int> num_bad_edges(3, 0);
        for (int i = 0; i < N_pts; ++i) {
          int i1 = i;
          int i2 = 0;
          if (i < N_pts - 1) {
            i2= i+1;
          }
          vec3_t x1, x2;
          m_Grid->GetPoint(pts[i1], x1.data());
          m_Grid->GetPoint(pts[i2], x2.data());
          double L = (x1 - x2).abs();
          vtkIdType id_neigh_cell = m_Part.c2cGG(id_cell, i);
          if (id_neigh_cell != -1) {
            bool bad_edge = false;
            vec3_t n2 = cellNormal(m_Grid, id_neigh_cell);
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
  cout << "done." << endl;
}

void FixCadGeometry::copyFaces(const QVector<bool> &copy_face)
{
  int N_del = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (!copy_face[id_cell]) {
      ++N_del;
    }
  }
  int N_cells = m_Grid->GetNumberOfCells() - N_del;
  int N_nodes = 0;
  QVector<bool> copy_node(m_Grid->GetNumberOfPoints(), false);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      if (copy_face[m_Part.n2cGG(id_node, i)]) {
        copy_node[id_node] = true;
        ++N_nodes;
        break;
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, N_cells, N_nodes);
  QVector<vtkIdType> old2new(m_Grid->GetNumberOfPoints(), -1);
  vtkIdType id_new_node = 0;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (copy_node[id_node]) {
      old2new[id_node] = id_new_node;
      vec3_t x;
      m_Grid->GetPoint(id_node, x.data());
      new_grid->GetPoints()->SetPoint(id_new_node, x.data());
      copyNodeData(m_Grid, id_node, new_grid, id_new_node);
      ++id_new_node;
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (copy_face[id_cell]) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      vtkIdType id_new_cell = new_grid->InsertNextCell(m_Grid->GetCellType(id_cell), N_pts, pts);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    }
  }
  makeCopy(new_grid, m_Grid);
  cout << "  deleted " << N_del << " faces" << endl;
}

void FixCadGeometry::fixNonManifold1()
{
  cout << "fixing non-manifold edges\n  (pass 1) ..." << endl;
  QVector<bool> copy_face(m_Grid->GetNumberOfCells(), true);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      if (m_Part.c2cGG(id_cell, i) == -1) {
        copy_face[id_cell] = false;
        break;
      }
    }
  }
  copyFaces(copy_face);
  cout << "done." << endl;
}

void FixCadGeometry::fixNonManifold2()
{
  cout << "fixing non-manifold edges\n  (pass 2) ..." << endl;
  QVector<bool> copy_face(m_Grid->GetNumberOfCells(), true);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    if (N_pts < 3) {
      copy_face[id_cell] = false;
      break;
    }
    QVector<QSet<vtkIdType> > n2c(N_pts);
    for (int i = 0; i < N_pts; ++i) {
      for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
        n2c[i].insert(m_Part.n2cGG(pts[i], j));
      }
    }
    QSet<vtkIdType> faces = n2c[0];
    for (int i = 1; i < N_pts; ++i) {
      faces = faces.intersect(n2c[i]);
    }
    if (faces.size() > 1) {
      vtkIdType id_del = id_cell;
      foreach (vtkIdType id, faces) {
        id_del = min(id_del, id);
      }
      copy_face[id_del] = false;
    }
  }
  copyFaces(copy_face);
  cout << "done." << endl;
}

void FixCadGeometry::operate()
{
  {
    DeleteStrayNodes del;
    del();
  }

  setAllCells();

  //prepare BCmap
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  QSet <int> bcs;
  l2g_t cells = m_Part.getCells();
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      bcs.insert(bc->GetValue(id_cell));
    }
  }
  QMap <int,int> BCmap;
  foreach(int bc, bcs) {
    BCmap[bc] = 1;
  }

  //set density infinite
  VertexMeshDensity VMD;
  VMD.density = 1e99;
  VMD.BCmap = BCmap;
  qWarning()<<"VMD.BCmap="<<VMD.BCmap;
  m_VMDvector.push_back(VMD);

  customUpdateNodeInfo();

  // fix non manifold edges
  fixNonManifold1();
  fixNonManifold2();

  //call surface mesher
  setGrid(m_Grid);
  setBoundaryCodes(bcs);
  setVertexMeshDensityVector(m_VMDvector);
  setDesiredLength();

  /*
  callMesher();

  // correct surface orientation
  CorrectSurfaceOrientation correct;
  correct.setGrid(m_Grid);
  correct.setAllCells();
  correct();
  */

  // finalise
  createIndices(m_Grid);
  customUpdateNodeInfo();
  setDesiredLength();
}

