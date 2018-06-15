// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "fixcadgeometry.h"
#include "surfacemesher.h"
#include "vertexmeshdensity.h"
#include "guimainwindow.h"
#include "deletestraynodes.h"
#include "correctsurfaceorientation.h"
#include "swaptriangles.h"
#include "insertpoints3d.h"
#include "deletecells.h"

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
  m_OriginalFeatureAngle = m_FeatureAngle;
  m_FeatureAngle = GeometryTools::deg2rad(0.5);
  m_AllowFeatureEdgeSwapping = false;
  m_AllowSmallAreaSwapping = true;
  m_NumDelaunaySweeps = 1;
  m_SnapTolerance = 1e-3;
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
  int count = 0;
  while (!done) {
    SwapTriangles swap;
    swap.setGrid(m_Grid);
    swap.setRespectBC(true);
    swap.setFeatureSwap(m_AllowFeatureEdgeSwapping);
    swap.setFeatureAngle(m_FeatureAngle);
    swap.setMaxNumLoops(1);
    swap.setSmallAreaSwap(m_AllowSmallAreaSwapping);
    swap.setDelaunayThreshold(1e6);
    swap.setVerboseOn();
    QSet<int> rest_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
    rest_bcs -= m_BoundaryCodes;
    swap.setBoundaryCodes(rest_bcs);
    swap();
    count = 0;
    if (swap.getNumSwaps() == 0 || count >= 100) {
      done = true;
    }
  }

  /*
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
  */
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
    node_type->SetValue(id_node, EG_FIXED_VERTEX);
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
            node_type->SetValue(pts[i], EG_SIMPLE_VERTEX);
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

void FixCadGeometry::markNonManifold()
{
  cout << "marking non-manifold edges\n" << endl;
  int new_bc = 0;
  m_NumNonManifold = 0;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QVector<bool> nm_face(m_Grid->GetNumberOfCells(), false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      new_bc = max(new_bc, cell_code->GetValue(id_cell));
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        QSet <vtkIdType> edge_cells;
        int N = 0;
        if (i < N_pts - 1) {
          N = getEdgeCells(pts[i], pts[i+1], edge_cells);
        } else {
          N = getEdgeCells(pts[i], pts[0], edge_cells);
        }
        if (N != 2) {
          nm_face[id_cell] = true;
          ++m_NumNonManifold;
          break;
        }
      }
    }
  }
  m_BoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();
  foreach (int bc, m_BoundaryCodes) {
    ++new_bc;
    QString bc_name = "NM_" + GuiMainWindow::pointer()->getBC(bc).getName();
    for (int level = 0; level < 2; ++level) {
      QVector<bool> new_nm_face = nm_face;
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        if (isSurface(id_cell, m_Grid)) {
          if ((cell_code->GetValue(id_cell)) == bc && nm_face[id_cell]) {
            vtkIdType N_pts, *pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            for (int i = 0; i < N_pts; ++i) {
              for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
                vtkIdType id_neigh = m_Part.n2cGG(pts[i], j);
                if (cell_code->GetValue(id_neigh) == bc ) {
                  new_nm_face[id_neigh] = true;
                }
              }
            }
          }
        }
      }
      nm_face = new_nm_face;
    }
    GuiMainWindow::pointer()->setBC(new_bc, BoundaryCondition(bc_name, "patch", new_bc));
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (cell_code->GetValue(id_cell) == bc && nm_face[id_cell]) {
          cell_code->SetValue(id_cell, new_bc);
        }
      }
    }
  }
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
  cout << "done." << endl;
}

void FixCadGeometry::createBox()
{
  vec3_t x1(1e20, 1e20, 1e20);
  vec3_t x2 = -1*x1;
  //
  EG_FORALL_NODES (id_node, m_Grid)
  {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      x1[i] = min(x[i], x1[i]);
      x2[i] = max(x[i], x2[i]);
    }
  }
  //
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  int num_new_cells = m_Grid->GetNumberOfCells() + 12;
  int num_new_nodes = m_Grid->GetNumberOfPoints() + 9;
  allocateGrid(new_grid, num_new_cells, num_new_nodes);
  makeCopyNoAlloc(m_Grid, new_grid);
  {
    vec3_t x[9];
    vec3_t dx = x2 - x1;
    //
    x[8] = 0.5*(x1 + x2);
    x1   = x[8] - dx;
    x2   = x[8] + dx;
    x[0] = vec3_t(x1[0], x1[1], x1[2]);
    x[1] = vec3_t(x2[0], x1[1], x1[2]);
    x[2] = vec3_t(x2[0], x2[1], x1[2]);
    x[3] = vec3_t(x1[0], x2[1], x1[1]);
    x[4] = vec3_t(x1[0], x1[1], x2[2]);
    x[5] = vec3_t(x2[0], x1[1], x2[2]);
    x[6] = vec3_t(x2[0], x2[1], x2[2]);
    x[7] = vec3_t(x1[0], x2[1], x2[1]);
    //
    EG_VTKDCN(vtkCharArray, node_type, new_grid, "node_type");
    for (int i = 0; i < 9; ++i) {
      vtkIdType id_node = m_Grid->GetNumberOfPoints() + i;
      new_grid->GetPoints()->SetPoint(id_node, x[i].data());
      if (i < 8) {
        node_type->SetValue(id_node, EG_OUTSIDE_VERTEX);
      } else {
        node_type->SetValue(id_node, EG_UNKNOWN_VERTEX);
      }
    }
    //
    QVector<QVector<vtkIdType> > pts(12);
    pts[0]  << 0 << 1 << 3 << 8;
    pts[1]  << 1 << 2 << 3 << 8;
    pts[2]  << 1 << 5 << 2 << 8;
    pts[3]  << 2 << 5 << 6 << 8;
    pts[4]  << 3 << 2 << 6 << 8;
    pts[5]  << 3 << 6 << 7 << 8;
    pts[6]  << 0 << 3 << 4 << 8;
    pts[7]  << 3 << 7 << 4 << 8;
    pts[8]  << 4 << 7 << 5 << 8;
    pts[9]  << 7 << 6 << 5 << 8;
    pts[10] << 0 << 4 << 1 << 8;
    pts[11] << 1 << 4 << 5 << 8;
    for (int i = 0; i < 12; ++i) {
      for (int j = 0; j < 4; ++j) {
        pts[i][j] += m_Grid->GetNumberOfPoints();
      }
      new_grid->InsertNextCell(VTK_TETRA, 4, pts[i].data());
    }
  }
  //
  makeCopy(new_grid, m_Grid);
}

void FixCadGeometry::computeCharLength()
{
  setAllCells();
  EG_VTKDCN(vtkDoubleArray, length, m_Grid, "node_meshdensity_desired");
  EG_FORALL_NODES (id_node, m_Grid) {
    double L = 1e99;
    vec3_t x1;
    m_Grid->GetPoint(id_node, x1.data());
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vec3_t x2;
      m_Grid->GetPoint(m_Part.n2nGG(id_node, i), x2.data());
      L = std::min(L, (x2 - x1).abs());
    }
    length->SetValue(id_node, L);
  }
}

FixCadGeometry::cut_t FixCadGeometry::snapCut(vtkIdType id_node1, vtkIdType id_node2)
{
  vec3_t x1, x2, n;
  cut_t C;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  //
  // determine if the edge cuts through the geometry
  //
  C.edge_cut = false;
  int count = 0;
  vec3_t x_old = 0.5*(x1 + x2);
  vec3_t x_snap;
  //
  while (!C.edge_cut) {
    //
    // snap and intersect
    //
    x_snap = m_Cad->snap(x_old, false);
    n      = m_Cad->getLastNormal();
    C.w    = intersection(x1, x2 - x1, x_snap, n);
    C.x    = x1 + C.w*(x2 - x1);
    //
    if ((C.w > 0 && C.w < 1) && ((C.x - x1).abs() >= m_SnapTolerance && (C.x - x2).abs() >= m_SnapTolerance)) {
      //
      // This edge potentially cuts through the surface
      //
      if ((x_old - C.x).abs() < m_SnapTolerance) {
        C.edge_cut = true;
      }
      ++count;
      x_old = C.x;
      //
      if (count > 10) {
        break;
      }
    } else {
      break;
    }
  }
  //
  C.R = m_Cad->getLastRadius();
  C.node1_surf = ( (x_snap - x1).abs() < m_SnapTolerance );
  C.node2_surf = ( (x_snap - x2).abs() < m_SnapTolerance );
  //
  vtkIdType id_face = m_Cad->getLastFaceId();
  EG_VTKDCN(vtkDoubleArray, length, m_Grid, "node_meshdensity_desired");
  EG_GET_CELL(id_face, m_Grid);
  C.L = length->GetValue(pts[0]);
  for (int i = 1; i < num_pts; ++i) {
    C.L = std::min(C.L, length->GetValue(pts[i]));
  }
  if (C.x[0] > 1.1 && C.edge_cut) {
    cout << "Baehh!!" << endl;
  }
  //
  return C;
}

void FixCadGeometry::marchOutside()
{
  setAllCells();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  //
  // mark outside nodes
  //
  bool changed = true;
  cout << "marching in ..." << endl;
  while (changed) {
    int N_surf = 0;
    int N_out  = 0;
    changed = false;
    EG_FORALL_NODES (id_node1, m_Grid) {
      if (node_type->GetValue(id_node1) == EG_OUTSIDE_VERTEX) {
        for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
          vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
          if (node_type->GetValue(id_node2) == EG_UNKNOWN_VERTEX) {
            cut_t C = snapCut(id_node1, id_node2);
            if (C.node1_surf) {
              EG_BUG;
            } else if (C.node2_surf) {
              node_type->SetValue(id_node2, EG_SURFACE_VERTEX);
              ++N_surf;
              changed = true;
            } else {
              if (!C.edge_cut) {
                node_type->SetValue(id_node2, EG_OUTSIDE_VERTEX);
                ++N_out;
                changed = true;
              }
            }
          }
        }
      }
    }
    cout << "  " << N_surf << " new surface nodes and " << N_out << " new outside nodes" << endl;
  }
  return;
  //
  // delete outside cells
  //
  DeleteCells del;
  del.setGrid(m_Grid);
  del.setAllVolumeCells();
  QList<vtkIdType> del_cells;
  EG_FORALL_CELLS (id_cell, m_Grid) {
    if (isVolume(id_cell, m_Grid)) {
      EG_GET_CELL_AND_TYPE(id_cell, m_Grid);
      if (type_cell != VTK_TETRA) {
        EG_BUG;
      }
      bool outside_node = false;
      for (int i = 0; i < num_pts; ++i) {
        if (node_type->GetValue(pts[i]) == EG_OUTSIDE_VERTEX) {
          outside_node = true;
          break;
        }
      }
      if (outside_node) {
        bool edge_cut = false;
        if (!edge_cut) edge_cut = snapCut(pts[0], pts[1]).edge_cut;
        if (!edge_cut) edge_cut = snapCut(pts[0], pts[2]).edge_cut;
        if (!edge_cut) edge_cut = snapCut(pts[0], pts[3]).edge_cut;
        if (!edge_cut) edge_cut = snapCut(pts[1], pts[2]).edge_cut;
        if (!edge_cut) edge_cut = snapCut(pts[1], pts[3]).edge_cut;
        if (!edge_cut) edge_cut = snapCut(pts[2], pts[3]).edge_cut;
        //
        if (!edge_cut) {
          del_cells << id_cell;
        }
      }
    }
  }
  del.setCellsToDelete(del_cells);
  del();
}

void FixCadGeometry::cut()
{
  setAllCells();
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  //
  // mark edges to cut
  //
  InsertPoints3d ins;
  ins.setGrid(m_Grid);
  ins.setAllVolumeCells();
  EG_FORALL_NODES (id_node1, m_Grid) {
    if (node_type->GetValue(id_node1) == EG_OUTSIDE_VERTEX) {
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        cut_t C = snapCut(id_node1, id_node2);
        if (C.edge_cut) {
          ins.addEdge(id_node1, id_node2, EG_SURFACE_VERTEX, C.x);
        }
      }
    }
  }
  ins();
}

void FixCadGeometry::refine()
{

}

void FixCadGeometry::operate()
{
  {
    DeleteStrayNodes del;
    del();
  }
  //
  GuiMainWindow::pointer()->storeCadInterfaces();
  m_Cad = dynamic_cast<CgalTriCadInterface*>(GuiMainWindow::pointer()->getCadInterface(0));
  if (!m_Cad) {
    EG_BUG;
  }
  //
  setAllCells();
  computeCharLength();
  createBox();

  for (int i = 0; i < 1; ++i) {
    marchOutside();
    cut();
  }
  marchOutside();
  updateCellIndex(m_Grid);
  updateNodeIndex(m_Grid);


  /*
  setAllCells();
  InsertPoints3d ins;
  ins.setGrid(m_Grid);
  ins.setAllVolumeCells();
  //
  for (int i = 0; i < m_Part.getNumberOfNodes(); ++i) {
    vtkIdType id_node1 = m_Part.globalNode(i);
    for (int j = 0; j < m_Part.n2nLSize(i); ++j) {
      vtkIdType id_node2 = m_Part.n2nLG(i, j);
      ins.addEdge(id_node1, id_node2, 0);
    }
  }
  //
  ins();
  */
}

