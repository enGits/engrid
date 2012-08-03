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
#include "booleangeometryoperation.h"
//#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include "guimainwindow.h"

#include <vtkGeometryFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkDelaunay3D.h>

void BooleanGeometryOperation::deleteNodes()
{
  QList<vtkIdType> keep_nodes;
  for (int i = 0; i < m_Part1.getNumberOfNodes(); ++i) {
    vtkIdType id_node = m_Part1.globalNode(i);
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    vec3_t x_proj = m_Proj2.project(x, -1, false);
    vec3_t n = m_Proj2.lastProjNormal();
    n.normalise();
    double L = 0;
    for (int j = 0; j < m_Part1.n2nGSize(id_node); ++j) {
      vec3_t xn;
      m_Grid->GetPoint(m_Part1.n2nGG(id_node, j), xn.data());
      L = max(L, (x-xn).abs());
    }
    if (m_Side2*((x - x_proj)*n) > 0.5*L) {
      keep_nodes.append(id_node);
    }
  }
  for (int i = 0; i < m_Part2.getNumberOfNodes(); ++i) {
    vtkIdType id_node = m_Part2.globalNode(i);
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    vec3_t x_proj = m_Proj1.project(x, -1, false);
    vec3_t n = m_Proj1.lastProjNormal();
    n.normalise();
    double L = 0;
    for (int j = 0; j < m_Part2.n2nGSize(id_node); ++j) {
      vec3_t xn;
      m_Grid->GetPoint(m_Part2.n2nGG(id_node, j), xn.data());
      L = max(L, (x-xn).abs());
    }
    if (m_Side1*((x - x_proj)*n) > 0.5*L) {
      keep_nodes.append(id_node);
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  makeCopy(m_Grid, grid);
  MeshPartition part(grid);
  part.setNodes(keep_nodes);
  part.extractToVtkGrid(m_Grid);
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
}

void BooleanGeometryOperation::checkOrientation()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  EG_FORALL_CELLS (id_face, m_Grid) {
    bool flip = false;
    if (m_Side1*m_Side2 < 0) {
      if (m_BCs1.contains(bc->GetValue(id_face)) && m_Side1 > 0) {
        flip = true;
      }
      if (m_BCs2.contains(bc->GetValue(id_face)) && m_Side2 > 0) {
        flip = true;
      }
    }
    if (flip) {
      EG_GET_CELL (id_face, m_Grid);
      QVector<vtkIdType> old_pts(num_pts);
      for (int i = 0; i < num_pts; ++i) {
        old_pts[i] = pts[i];
      }
      for (int i = 0; i < num_pts; ++i) {
        pts[i] = old_pts[num_pts - i - 1];
      }
    }
  }
}

bool BooleanGeometryOperation::fillGap_prepare()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  bool found = false;
  int N = m_Grid->GetNumberOfPoints();
  m_OpenNode.fill(0, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      QVector<vtkIdType> ids(N_pts + 1);
      for (int i = 0; i <= N_pts; ++i) {
        ids[i] = pts[i];
      }
      ids[N_pts] = ids[0];
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        if (m_Part.c2cGG(id_cell, i) == -1) {
          int on = 0;
          QList<int> bcs1 = m_BCs1;
          QList<int> bcs2 = m_BCs2;
          int bc = cell_code->GetValue(id_cell);
          if (m_BCs1.contains(cell_code->GetValue(id_cell))) {
            on = 1;
          } else if (m_BCs2.contains(cell_code->GetValue(id_cell))) {
            on = 2;
          } else {
            EG_BUG;
          }
          m_OpenNode[ids[i]]   = on;
          m_OpenNode[ids[i+1]] = on;
          if (!found) {
            m_CurrentEdge.id1 = ids[i+1];
            m_CurrentEdge.id2 = ids[i];
            m_CurrentEdge.bc = cell_code->GetValue(id_cell);
            m_CurrentTriangle.bc = m_CurrentEdge.bc;
            m_CurrentTriangle.id1 = ids[0];
            m_CurrentTriangle.id2 = ids[1];
            m_CurrentTriangle.id3 = ids[2];
            found = true;
          }
        }
      }
    }
  }
  if (!found) {
    return false;
  }
  m_Node2Cell.clear();
  m_Node2Cell.resize(m_Grid->GetNumberOfPoints());
  m_Node2Node.clear();
  m_Node2Node.resize(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_OpenNode[id_node]) {
      m_SplitNode[id_node] = true;
      m_Node2Cell[id_node].clear();
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        m_Node2Cell[id_node].insert(m_Part.n2cGG(id_node, i));
      }
      m_Node2Node[id_node].clear();
      for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
        m_Node2Node[id_node].insert(m_Part.n2nGG(id_node, i));
      }
    }
  }
  return true;
}

void BooleanGeometryOperation::fillGap_updateOpenNode(vtkIdType id_node)
{
  if (m_OpenNode[id_node]) {
    if (id_node == 571) {
      cout << "break-point" << endl;
    }
    int on = 0;
    foreach (vtkIdType id_neigh, m_Node2Node[id_node]) {
      QSet<vtkIdType> shared_cells = m_Node2Cell[id_node];
      shared_cells.intersect(m_Node2Cell[id_neigh]);
      if (shared_cells.size() == 1) {
        on = m_OpenNode[id_node];
        break;
      }
    }
    m_OpenNode[id_node] = on;
  }
}

bool BooleanGeometryOperation::fillGap_onDifferentSides(vtkIdType id1, vtkIdType id2, vtkIdType id3)
{
  if (m_OpenNode[id1] && m_OpenNode[id2] && m_OpenNode[id3]) {
    if (m_OpenNode[id1] != m_OpenNode[id2]) {
      return true;
    }
    if (m_OpenNode[id1] != m_OpenNode[id3]) {
      return true;
    }
    if (m_OpenNode[id2] != m_OpenNode[id3]) {
      return true;
    }
  }
  return false;
}

bool BooleanGeometryOperation::fillGap_step()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QSet<vtkIdType> candidates;
  vtkIdType id1 = m_CurrentEdge.id1;
  vtkIdType id2 = m_CurrentEdge.id2;
  if (m_Triangles.size() == 0) {
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_OpenNode[id_node]) {
        QSet<int> bcs;
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          bcs.insert(cell_code->GetValue(m_Part.n2cGG(id_node, i)));
        }
        if (!bcs.contains(m_CurrentEdge.bc)) {
          candidates.insert(id_node);
        }
      }
    }
  } else {
    foreach (vtkIdType id_node, m_Node2Node[id1]) {
      if (m_OpenNode[id_node] == m_OpenNode[id1]) {
        candidates.insert(id_node);
      }
    }
    foreach (vtkIdType id_node, m_Node2Node[id2]) {
      if (m_OpenNode[id_node] == m_OpenNode[id2]) {
        candidates.insert(id_node);
      }
    }
  }
  if (candidates.size() == 0) {
    return false;
  }
  vec3_t n1 = triNormal(m_Grid, m_CurrentTriangle.id1, m_CurrentTriangle.id2, m_CurrentTriangle.id3);
  double A1 = n1.abs();
  double d_min = 1e99;
  double amm = 1e99; // maximal minimal angle
  vec3_t x1, x2;
  m_Grid->GetPoint(id1, x1.data());
  m_Grid->GetPoint(id2, x2.data());
  tri_t T;
  bool found = false;
  foreach (vtkIdType id3, candidates) {
    vec3_t n2 = triNormal(m_Grid, id1, id2, id3);
    double A2 = n2.abs();
    if (n1*n2 > 0 || m_Triangles.size() == 0) {
      vec3_t x3;
      m_Grid->GetPoint(id3, x3.data());
      double d = (x3 - 0.5*(x1+x2)).abs();
      vec3_t v1 = x2 - x1;
      vec3_t v2 = x3 - x2;
      vec3_t v3 = x1 - x3;
      double a1 = angle(v1, v2);
      double a2 = angle(v2, v3);
      double a3 = angle(v3, v1);
      double a_min = min(a1, min(a2, a3));
      if (d < d_min) {
        d_min = d;
        //amm = a_min;
        T.id1 = id1;
        T.id2 = id2;
        T.id3 = id3;
        T.bc = m_CurrentEdge.bc;
        found = true;
      }
    }
  }
  if (found) {
    m_Triangles.append(T);
  } else {
    return false;
  }
  m_CurrentTriangle = T;

  // update node -> cell structure
  m_Node2Cell[T.id1].insert(-m_Triangles.size());
  m_Node2Cell[T.id2].insert(-m_Triangles.size());
  m_Node2Cell[T.id3].insert(-m_Triangles.size());

  // update node -> node structure
  m_Node2Node[T.id1].insert(T.id2);
  m_Node2Node[T.id1].insert(T.id3);
  m_Node2Node[T.id2].insert(T.id1);
  m_Node2Node[T.id2].insert(T.id3);
  m_Node2Node[T.id3].insert(T.id1);
  m_Node2Node[T.id3].insert(T.id2);

  // update open node structure
  int open1 = m_OpenNode[T.id1];
  int open2 = m_OpenNode[T.id2];
  int open3 = m_OpenNode[T.id3];
  fillGap_updateOpenNode(T.id1);
  fillGap_updateOpenNode(T.id2);
  fillGap_updateOpenNode(T.id3);

  int is_open1 = m_OpenNode[T.id1];
  int is_open2 = m_OpenNode[T.id2];
  int is_open3 = m_OpenNode[T.id3];
  if        (is_open2 && is_open3) {
    m_CurrentEdge.id1 = T.id3;
    m_CurrentEdge.id2 = T.id2;
  } else if (is_open3 && is_open1) {
    m_CurrentEdge.id1 = T.id1;
    m_CurrentEdge.id2 = T.id3;
  } else {
    return false;
  }

  return true;
}

void BooleanGeometryOperation::fillGap_createTriangles()
{
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, m_Grid->GetNumberOfCells() + m_Triangles.size(), m_Grid->GetNumberOfPoints());
  makeCopyNoAlloc(m_Grid, grid);
  {
    EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
    vtkIdType pts[3];
    foreach (tri_t T, m_Triangles) {
      pts[0] = T.id1;
      pts[1] = T.id2;
      pts[2] = T.id3;
      vtkIdType id_new = grid->InsertNextCell(VTK_TRIANGLE, 3, pts);
      bc->SetValue(id_new, T.bc);
    }
  }
  makeCopy(grid, m_Grid);
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
}

void BooleanGeometryOperation::fillGap()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  m_SplitNode.fill(false, m_Grid->GetNumberOfPoints());
  do {
    m_Triangles.clear();
    if (fillGap_prepare()) {
      int count = 0;
      bool done = false;
      while (count < 1000 && !done) {
        if (!fillGap_step()) {
          done = true;
        }
        done = true;
        ++count;
      }
      fillGap_createTriangles();
    }
    m_Triangles.clear();
  } while (m_Triangles.size() > 0);
}

void BooleanGeometryOperation::smoothJunction_triangulate()
{
  QVector<bool> split_cell(m_Grid->GetNumberOfCells(), false);
  for (int layer = 0; layer < m_NumCutLayers; ++layer) {
    EG_FORALL_NODES(id_node, m_Grid) {
      if (m_SplitNode[id_node]) {
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          vtkIdType id_cell = m_Part.n2cGG(id_node, i);
          if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
            split_cell[m_Part.n2cGG(id_node, i)] = true;
          }
        }
      }
    }
    EG_FORALL_CELLS(id_cell, m_Grid) {
      if (split_cell[id_cell]) {
        EG_GET_CELL(id_cell, m_Grid);
        for (int i = 0; i < num_pts; ++i) {
          m_SplitNode[pts[i]] = true;
        }
      }
    }
  }
  int num_new_faces = 0;
  foreach (bool split, split_cell) {
    if (split) {
      ++num_new_faces;
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, m_Grid->GetNumberOfCells() + num_new_faces, m_Grid->GetNumberOfPoints());
  EG_FORALL_NODES(id_node, m_Grid) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    grid->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(m_Grid, id_node, grid, id_node);
  }
  EG_FORALL_CELLS(id_cell, m_Grid) {
    EG_GET_CELL(id_cell, m_Grid) {
      if (split_cell[id_cell]) {
        vtkIdType tri_pts[3];
        tri_pts[0] = pts[0];
        tri_pts[1] = pts[1];
        tri_pts[2] = pts[2];
        vtkIdType id_new;
        id_new = grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
        copyCellData(m_Grid, id_cell, grid, id_new);
        tri_pts[0] = pts[2];
        tri_pts[1] = pts[3];
        tri_pts[2] = pts[0];
        id_new = grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
        copyCellData(m_Grid, id_cell, grid, id_new);
      } else {
        vtkIdType id_new = grid->InsertNextCell(type_cell, num_pts, pts);
        copyCellData(m_Grid, id_cell, grid, id_new);
      }
    }
  }
  makeCopy(grid, m_Grid);
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
}

void BooleanGeometryOperation::smoothJunction_updateBCs()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QList<int> bcs = m_BCs1 + m_BCs2;
  EG_FORALL_CELLS(id_cell, m_Grid) {
    EG_GET_CELL(id_cell, m_Grid);
    if (type_cell == VTK_TRIANGLE) {
      int bc_best = -1;
      double L_min = 1e99;
      foreach (int bc, bcs) {
        double L_max = 0;
        for (int i = 0; i < num_pts; ++i) {
          vec3_t x_node;
          m_Grid->GetPoint(pts[i], x_node.data());
          vec3_t x_proj = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_node);
          L_max = max(L_max, (x_node - x_proj).abs());
        }
        if (L_max < L_min) {
          bc_best = bc;
          L_min = L_max;
        }
      }
      if (bc_best == -1) {
        EG_BUG;
      }
      cell_code->SetValue(id_cell, bc_best);
    }
  }
}

double BooleanGeometryOperation::smoothJunction_mesher()
{
  prepare();
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
  updateNodeInfo(true);
  computeMeshDensity();
  int inserted_nodes = insertNodes();
  int num_nodes = m_Grid->GetNumberOfPoints();
  updateNodeInfo();
  swap();
  int deleted_nodes = SurfaceAlgorithm::deleteNodes();
  for (int i = 0; i < 2; ++i) {
    smooth(1, true);
    swap();
  }
  createIndices(m_Grid);
  updateNodeInfo(false);
  return double(inserted_nodes - deleted_nodes)/double(num_nodes);
}

void BooleanGeometryOperation::smoothJunction()
{
  smoothJunction_triangulate();
  smoothJunction_updateBCs();
  double change_ratio;
  for (int iter = 1; iter <= 20; ++iter) {
    change_ratio = smoothJunction_mesher();
    smoothJunction_updateBCs();
    cout << "iteration " << iter << ": change ratio = " << change_ratio*100 << "%" << endl;
    if (change_ratio < 0.01 && iter >= 5) {
      break;
    }
  }
}

void BooleanGeometryOperation::operate()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  deleteNodes();
  checkOrientation();
  fillGap();
  //smoothJunction();
}
