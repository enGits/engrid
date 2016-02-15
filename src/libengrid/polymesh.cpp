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

#include "polymesh.h"
#include "polymolecule.h"
#include "guimainwindow.h"

#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCurvatures.h>

// ==========
//   face_t
// ==========

PolyMesh::face_t::face_t(int N, int o, int n, vec3_t rv, int b)
{
  node.resize(N);
  owner = o;
  neighbour = n;
  bc = b;
  ref_vec = rv;
}


// ==========
//   node_t
// ==========

PolyMesh::node_t::node_t(const QVector<vtkIdType> &ids)
{
  id = ids;
  qSort(id);
}

PolyMesh::node_t::node_t(vtkIdType id1)
{
  id.resize(1);
  id[0] = id1;
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2)
{
  id.resize(2);
  id[0] = id1;
  id[1] = id2;
  qSort(id);
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3)
{
  id.resize(3);
  id[0] = id1;
  id[1] = id2;
  id[2] = id3;
  qSort(id);
}

PolyMesh::node_t::node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType id4)
{
  id.resize(4);
  id[0] = id1;
  id[1] = id2;
  id[2] = id3;
  id[3] = id4;
  qSort(id);
}

bool PolyMesh::node_t::operator<(const PolyMesh::node_t &N) const
{
  int num = min(id.size(), N.id.size());
  for (int i = 0; i < num; ++i) {
    if (id[i] < N.id[i]) {
      return true;
    }
    if (id[i] > N.id[i]) {
      return false;
    }
  }
  if (id.size() < N.id.size()) {
    return true;
  }
  return false;
}

bool PolyMesh::node_t::operator>(const PolyMesh::node_t &N) const
{
  int num = min(id.size(), N.id.size());
  for (int i = 0; i < num; ++i) {
    if (id[i] > N.id[i]) {
      return true;
    }
    if (id[i] < N.id[i]) {
      return false;
    }
  }
  if (id.size() > N.id.size()) {
    return true;
  }
  return false;
}

bool PolyMesh::node_t::operator==(const PolyMesh::node_t &N) const
{
  if (id.size() != N.id.size()) {
    return false;
  }
  for (int i = 0; i < id.size(); ++i) {
    if (id[i] != N.id[i]) {
      return false;
    }
  }
  return true;
}





// ============
//   PolyMesh
// ============

PolyMesh::PolyMesh(vtkUnstructuredGrid *grid, bool dualise, double pull_in, bool optimise, bool split_faces, bool split_cells)
{
  m_CreateDualMesh = false;
  m_OptimiseConvexity = false;
  m_SplitCells = false;
  m_SplitFaces = split_faces;
  if (dualise) {
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      if (isVolume(id_cell, grid) && grid->GetCellType(id_cell) != VTK_POLYHEDRON) {
        m_CreateDualMesh = true;
        break;
      }
    }
    if (m_CreateDualMesh) {
      m_OptimiseConvexity = optimise;

    }
  }

  m_AttractorWeight = 0.0;
  m_PullInFactor = pull_in;
  m_OptimiseConvexity = optimise && m_CreateDualMesh;
  m_Grid = grid;
  m_Part.setGrid(m_Grid);
  m_Part.setAllCells();
  findPolyCells();
  createNodesAndFaces();
  checkFaceOrientation();
  buildPoint2Face();
  buildPCell2Face();
  computePoints();

  if (m_OptimiseConvexity) {
    for (int iter = 0; iter < 2; ++iter) {
      int num_bad = 0;
      int i_improve = 0;
      for (int i = 0; i < numCells(); ++i) {

        // check if any of the faces is a boundary face
        bool no_boundary = true;
        for (int j = 0; j < numFacesOfPCell(i); ++j) {
          int bc = boundaryCode(pcell2Face(i, j));
          if (bc != 0) {
            no_boundary = false;
            break;
          }
        }

        if (no_boundary) {
          PolyMolecule pm(this, i);
          if (!pm.allPositive()) {
            ++i_improve;
            pm.optimise();
            if (pm.minPyramidVolume() < 0) {
              ++num_bad;
            }
          }
        }
      }
      cout << i_improve << " cells out of " << numCells() << " were optimised." << endl;
      if (num_bad == 0) {
        break;
      }
    }
  }

  if (m_SplitFaces) {
    splitConcaveFaces();
  }

  if (m_SplitCells) {
    for (int iter = 0; iter < 0; ++iter) {
      int num_improved = 0;
      int num_cells = numCells();
      for (int i = 0; i < num_cells; ++i) {
        PolyMolecule pm(this, i);
        if (!pm.allPositive()) {
          ++num_improved;
          pm.centreSplit();
        }
      }
      buildPoint2Face();
      buildPCell2Face();
      cout << num_improved << " cells out of " << num_cells << " were split." << endl;
      if (num_improved == 0) {
        break;
      }
    }
  }
  sortFaces();
  buildPoint2Face();
  buildPCell2Face();
}

void PolyMesh::merge(PolyMesh *poly)
{
  foreach (face_t face, poly->m_Faces) {
    for (int i = 0; i < face.node.size(); ++i) {
      face.node[i] += m_Points.size();
    }
    face.owner += m_NumPolyCells;
    face.neighbour += m_NumPolyCells;
    m_Faces.append(face);
  }
  m_Points += poly->m_Points;
  m_NumPolyCells += poly->m_NumPolyCells;
  sortFaces();
  collectBoundaryConditions();
}

void PolyMesh::triangulateBadFaces()
{
  m_IsBadCell.fill(false, numCells());
  int num_bad = 0;
  for (int i = 0; i < numCells(); ++i) {
    PolyMolecule pm(this, i);
    if (!pm.allPositive()) {
      ++num_bad;
      m_IsBadCell[i] = true;
    }
  }
  QList<int> bad_faces;
  for (int i = 0; i < numFaces(); ++i) {
    if (m_IsBadCell[owner(i)]) {
      bad_faces.append(i);
    } else if (neighbour(i) != -1) {
      if (m_IsBadCell[neighbour(i)]) {
        bad_faces.append(i);
      }
    }
  }

  foreach (int i_face, bad_faces) {
    face_t face = m_Faces[i_face];
    QVector<vec3_t> x(face.node.size());
    for (int i = 0; i < x.size(); ++i) {
      x[i] = nodeVector(face.node[i]);
    }
    EG_VTKSP(vtkPolyData, poly);
    createPolyData(x, poly);
    EG_VTKSP(vtkTriangleFilter, tri);
    tri->SetInputData(poly);
    tri->Update();
    if (tri->GetOutput()->GetNumberOfPoints() > face.node.size()) {
      EG_BUG;
    }
    QVector<face_t> new_faces(tri->GetOutput()->GetNumberOfCells(), face);
    for (int i = 0; i < new_faces.size(); ++i) {
      new_faces[i].node.resize(3);
      EG_GET_CELL(i, tri->GetOutput());
      for (int j = 0; j < 3; ++j) {
        new_faces[i].node[j] = face.node[pts[j]];
      }
    }
    m_Faces[i_face] = new_faces[0];
    for (int i = 1; i < new_faces.size(); ++i) {
      m_Faces.append(new_faces[i]);
    }
  }

  sortFaces();
  buildPoint2Face();
  buildPCell2Face();

  cout << num_bad << " cells out of " << numCells() << " are concave" << endl;
  cout << bad_faces.size() << " faces have been triangulated" << endl;
}

/*
void PolyMesh::splitConcaveFaces()
{
  for (int i_cell = 0; i_cell < numCells(); ++i_cell) {
    if (m_IsBadCell[i_cell]) {
      PolyMolecule pm(this, i_cell);
      EG_VTKSP(vtkPolyData, poly);
      pm.createPolyData(poly);
      EG_VTKSP(vtkCurvatures, curv);
      curv->SetCurvatureTypeToMinimum();
      curv->SetInput(poly);
      curv->Update();
      EG_VTKDCN(vtkDoubleArray, curvature, curv->GetOutput(), "Minimum_Curvature");
      vtkIdType id_ref = 0;
      vec3_t x_ref;
      double curv_min = 1e99;
      EG_FORALL_NODES(id_node, poly) {
        if (curvature->GetValue(id_node) < curv_min) {
          curv_min = curvature->GetValue(id_node);
          id_ref = id_node;
          poly->GetPoint(id_ref, x_ref.data());
        }
      }
      vec3_t n_ref;
      EG_FORALL_CELLS(id_face, poly) {
        EG_GET_CELL(id_face, poly);
        bool found = false;
        if (type_cell != VTK_TRIANGLE) {
          EG_BUG;
        }
        for (int i = 0; i < num_pts; ++i) {
          if (pts[i] == id_ref) {
            found = true;
            break;
          }
        }
        if (found) {
          QVector<vec3_t> x(num_pts);
          for (int i = 0; i < num_pts; ++i) {
            poly->GetPoint(pts[i], x[i].data());
          }
          n_ref = (x[1] - x[0]);
          n_ref = n_ref.cross(x[2] - x[0]);
          n_ref.normalise();
          break;
        }
      }


      EG_VTKSP(vtkDoubleArray, div);
      div->SetNumberOfComponents(1);
      div->SetNumberOfTuples(poly->GetNumberOfPoints());
      div->SetName("div");

      bool done = false;
      vtkIdType id_node0 = id_ref;
      vtkIdType id_node1 = id_ref;
      QVector<QSet<vtkIdType> > n2n;
      QVector<QSet<vtkIdType> > n2c;
      QVector<QVector<vtkIdType> > c2c;
      createPolyDataN2N(poly, n2n);
      createPolyDataN2C(poly, n2c);
      createPolyDataC2C(poly, c2c);

      vec3_t n_node(0,0,0);
      foreach (vtkIdType id_face, n2c[id_ref]) {
        EG_GET_CELL(id_face, poly);
        QVector<vec3_t> x(3);
        for (int i = 0; i < 3; ++i) {
          poly->GetPoint(pts[i], x[i].data());
        }
        n_node += triNormal(x[0], x[1], x[2]);
      }
      n_node.normalise();

      QVector<bool> div_node(poly->GetNumberOfPoints(), false);
      EG_FORALL_NODES(id_node, poly) {
        div->SetValue(id_node, 0);
      }
      div_node[id_ref] = true;
      double dv = 1.0;
      div->SetValue(id_ref, dv);
      dv += 0.2;
      int count = 0;
      double c_min = 1e99;
      double L = 1e99;
      foreach (vtkIdType id, n2n[id_node0]) {
        vec3_t x;
        poly->GetPoint(id, x.data());
        L = min(L, (x - x_ref).abs());
        double c = curvature->GetValue(id);
        if (c < c_min) {
          c_min = c;
          id_node1 = id;
        }
      }
      div_node[id_node1] = true;
      div->SetValue(id_node1, dv);
      dv += 0.2;
      while (!done) {
        vtkIdType id_node2 = id_node1;
        double d_min = 1e99;
        foreach (vtkIdType id, n2n[id_node1]) {
          if (id == id_ref && id_node0 != id_ref) {
            done = true;
            break;
          }
        }
        if (!done) {
          foreach (vtkIdType id, n2n[id_node1]) {
            if (!div_node[id]) {
              vec3_t x;
              poly->GetPoint(id, x.data());
              if ((x - x_ref)*n_node < 0.25*L) {
                double d = (x - x_ref)*n_ref;
                if (fabs(d) < d_min && d < 0.1*L) {
                  d_min = d;
                  id_node2 = id;
                }
              }
            }
          }
          div_node[id_node2] = true;
          div->SetValue(id_node2, dv);
          dv += 0.2;
          id_node0 = id_node1;
          id_node1 = id_node2;
          ++count;
          if (count == 100) {
            EG_BUG;
          }
        }
      }
      poly->GetPointData()->AddArray(div);
      EG_VTKSP(vtkXMLPolyDataWriter, vtp);
      vtp->SetFileName(qPrintable(GuiMainWindow::pointer()->getCwd() + "/triangulated_cell.vtp"));
      vtp->SetInput(poly);
      vtp->SetDataModeToAscii();
      vtp->Write();
      EG_BUG;
    }
  }
}
*/

void PolyMesh::getFacesOfEdgeInsideCell(vtkIdType id_cell, vtkIdType id_node1, vtkIdType id_node2, int &face1, int &face2)
{
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell, num_pts, pts);
  face1 = -1;
  face2 = -1;

  if (m_Grid->GetCellType(id_cell) == VTK_TETRA) {
    if        (id_node1 == pts[0]) {
      if      (id_node2 == pts[1]) { face1 = 0; face2 = 1; }
      else if (id_node2 == pts[2]) { face1 = 2; face2 = 0; }
      else if (id_node2 == pts[3]) { face1 = 1; face2 = 2; }
    } else if (id_node1 == pts[1]) {
      if      (id_node2 == pts[0]) { face1 = 1; face2 = 0; }
      else if (id_node2 == pts[2]) { face1 = 0; face2 = 3; }
      else if (id_node2 == pts[3]) { face1 = 3; face2 = 1; }
    } else if (id_node1 == pts[2]) {
      if      (id_node2 == pts[0]) { face1 = 0; face2 = 2; }
      else if (id_node2 == pts[1]) { face1 = 3; face2 = 0; }
      else if (id_node2 == pts[3]) { face1 = 2; face2 = 3; }
    } else if (id_node1 == pts[3]) {
      if      (id_node2 == pts[0]) { face1 = 2; face2 = 1; }
      else if (id_node2 == pts[1]) { face1 = 1; face2 = 3; }
      else if (id_node2 == pts[2]) { face1 = 3; face2 = 2; }
    }
  } else if (m_Grid->GetCellType(id_cell) == VTK_PYRAMID) {
    if        (id_node1 == pts[0]) {
      if      (id_node2 == pts[1]) { face1 = 0; face2 = 1; }
      else if (id_node2 == pts[3]) { face1 = 4; face2 = 0; }
      else if (id_node2 == pts[4]) { face1 = 1; face2 = 4; }
    } else if (id_node1 == pts[1]) {
      if      (id_node2 == pts[0]) { face1 = 1; face2 = 0; }
      else if (id_node2 == pts[2]) { face1 = 0; face2 = 2; }
      else if (id_node2 == pts[4]) { face1 = 2; face2 = 1; }
    } else if (id_node1 == pts[2]) {
      if      (id_node2 == pts[1]) { face1 = 2; face2 = 0; }
      else if (id_node2 == pts[3]) { face1 = 0; face2 = 3; }
      else if (id_node2 == pts[4]) { face1 = 3; face2 = 2; }
    } else if (id_node1 == pts[3]) {
      if      (id_node2 == pts[0]) { face1 = 0; face2 = 4; }
      else if (id_node2 == pts[2]) { face1 = 3; face2 = 0; }
      else if (id_node2 == pts[4]) { face1 = 4; face2 = 3; }
    } else if (id_node1 == pts[4]) {
      if      (id_node2 == pts[0]) { face1 = 4; face2 = 1; }
      else if (id_node2 == pts[1]) { face1 = 1; face2 = 2; }
      else if (id_node2 == pts[2]) { face1 = 2; face2 = 3; }
      else if (id_node2 == pts[3]) { face1 = 3; face2 = 4; }
    }
  } else {
    EG_BUG;
  }

  if (face1 == -1 || face2 == -1) {
    EG_BUG;
  }
}

void PolyMesh::getSortedEdgeCells(vtkIdType id_node1, vtkIdType id_node2, QList<vtkIdType> &cells, bool &is_loop)
{
  cells.clear();
  vtkIdType id_start = -1;
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, i);
    if (isVolume(id_cell, m_Grid)) {
      if (m_Cell2PCell[id_cell] == -1) {
        QVector<vtkIdType> pts;
        getPointsOfCell(m_Grid, id_cell, pts);
        if (pts.contains(id_node2)) {
          id_start = id_cell;
        }
      }
    }
  }
  if (id_start == -1) {
    EG_BUG;
  }
  vtkIdType id_cell_old = id_start;
  vtkIdType id_cell_new = id_start;
  cells.append(id_cell_new);
  bool finished = false;
  bool prepend = false;
  while (!finished) {
    do {
      int f1, f2;
      getFacesOfEdgeInsideCell(id_cell_new, id_node1, id_node2, f1, f2);
      vtkIdType id_neigh1 = m_Part.c2cGG(id_cell_new, f1);
      vtkIdType id_neigh2 = m_Part.c2cGG(id_cell_new, f2);
      if (id_neigh1 == -1) EG_BUG;
      if (id_neigh2 == -1) EG_BUG;
      if (id_neigh1 == id_cell_old) {
        id_cell_old = id_cell_new;
        id_cell_new = id_neigh2;
      } else {
        id_cell_old = id_cell_new;
        id_cell_new = id_neigh1;
      }
      if (prepend) {
        cells.prepend(id_cell_new);
      } else {
        cells.append(id_cell_new);
      }
    } while (id_cell_new != id_start  &&  !isSurface(id_cell_new, m_Grid)  &&  m_Cell2PCell[id_cell_new] == -1);
    if ((isSurface(id_cell_new, m_Grid) || m_Cell2PCell[id_cell_new] != -1) && !prepend) {
      id_cell_new = cells[0];
      id_cell_old = cells[1];
      prepend = true;
    } else {
      finished = true;
    }
  }

  // remove last cell for loops
  if (cells.size() > 1 && cells.first() == cells.last() && m_Cell2PCell[cells.first()] == -1) {
    cells.pop_back();
    is_loop = true;
  } else {
    is_loop = false;
  }
}

bool PolyMesh::isDualFace(vtkIdType id_face)
{
  vtkIdType id_cell = m_Part.getVolumeCell(id_face);
  if (m_Cell2PCell[id_cell] == -1) {
    return true;
  }
  return false;
}

void PolyMesh::getSortedPointFaces(vtkIdType id_node, int bc, QList<vtkIdType> &faces, bool &is_loop)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  faces.clear();
  vtkIdType id_start = -1;
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node, i);
    if (isSurface(id_cell, m_Grid)) {
      if (cell_code->GetValue(id_cell) == bc && isDualFace(id_cell)) {
        id_start = id_cell;
        break;
      }
    }
  }
  if (id_start == -1) {
    return;
  }
  vtkIdType id_face_old = id_start;
  vtkIdType id_face_new = id_start;
  faces.append(id_face_new);
  bool finished = false;
  bool prepend = false;
  while (!finished) {
    do {
      QList<vtkIdType> id_neigh;
      for (int i = 0; i < m_Part.c2cGSize(id_face_new); ++i) {
        if (cellContainsNode(m_Grid, m_Part.c2cGG(id_face_new, i), id_node)) {
          id_neigh.append(m_Part.c2cGG(id_face_new, i));
        }
      }
      if (id_neigh[0] == id_face_old) {
        id_face_old = id_face_new;
        id_face_new = id_neigh[1];
      } else {
        id_face_old = id_face_new;
        id_face_new = id_neigh[0];
      }
      if (cell_code->GetValue(id_face_new) == bc  &&  isDualFace(id_face_new)) {
        if (prepend) {
          faces.prepend(id_face_new);
        } else {
          faces.append(id_face_new);
        }
      }
    } while (id_face_new != id_start  &&  cell_code->GetValue(id_face_new) == bc  &&  isDualFace(id_face_new));
    if ((cell_code->GetValue(id_face_new) != bc || !isDualFace(id_face_new))  &&  !prepend) {
      id_face_old = id_face_new;
      id_face_new = faces[0];
      if (faces.size() > 1) {
        id_face_old = faces[1];
      }
      prepend = true;
    } else {
      finished = true;
    }
  }

  // remove last face for loops
  if (faces.size() > 1 && faces.first() == faces.last()) {
    faces.pop_back();
    is_loop = true;
  } else {
    is_loop = false;
  }

}

void PolyMesh::findPolyCells()
{
  m_Cell2PCell.fill(-1, m_Grid->GetNumberOfCells());
  m_Node2PCell.fill(-1, m_Grid->GetNumberOfPoints());
  m_NumPolyCells = 0;
  if (m_CreateDualMesh) {
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!isHexCoreNode(id_node)) {
        bool tetra_or_pyramid = false;
        for (int j = 0; j < m_Part.n2cGSize(id_node); ++j) {
          vtkIdType id_cell = m_Part.n2cGG(id_node, j);
          vtkIdType type_cell = m_Grid->GetCellType(id_cell);
          if (type_cell == VTK_TETRA || type_cell == VTK_PYRAMID) {
            tetra_or_pyramid = true;
          }
        }
        if (tetra_or_pyramid) {
          m_Node2PCell[id_node] = m_NumPolyCells;
          ++m_NumPolyCells;
        }
      }
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if (isVolume(id_cell, m_Grid)) {
      if (type_cell == VTK_WEDGE || type_cell == VTK_HEXAHEDRON || type_cell == VTK_POLYHEDRON || !m_CreateDualMesh) {
        m_Cell2PCell[id_cell] = m_NumPolyCells;
        ++m_NumPolyCells;
      }
    }
  }
  m_CellCentre.resize(m_NumPolyCells);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_Node2PCell[id_node] != -1) {
      m_Grid->GetPoint(id_node, m_CellCentre[m_Node2PCell[id_node]].data());
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Cell2PCell[id_cell] != -1) {
      m_CellCentre[m_Cell2PCell[id_cell]] = cellCentre(m_Grid, id_cell);
    }
  }
}

void PolyMesh::createFace(QList<node_t> nodes, int owner, int neighbour, vec3_t ref_vec, int bc)
{
  if (owner > neighbour && neighbour != -1) {
    swap(owner, neighbour);
    ref_vec *= -1;
  }
  face_t face(nodes.size(), owner, neighbour, ref_vec, bc);
  for (int i = 0; i < nodes.size(); ++i) {
    int idx = m_Nodes.insert(nodes[i]);
    face.node[i] = idx;
  }
  m_Faces.append(face);
}

void PolyMesh::createCornerFace(vtkIdType id_cell, int i_face, vtkIdType id_node)
{
  QList<vtkIdType> edge_nodes;
  QVector<vtkIdType> face_nodes;
  getFaceOfCell(m_Grid, id_cell, i_face, face_nodes);
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints(id_cell, num_pts, pts);
  for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
    vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
    if (face_nodes.contains(id_neigh)) {
      edge_nodes.append(id_neigh);
    }
  }
  if (edge_nodes.size() != 2) {
    EG_BUG;
  }
  int owner = m_Cell2PCell[id_cell];
  if (owner == -1) {
    EG_BUG;
  }
  int neighbour = m_Node2PCell[id_node];
  int bc = 0;
  if (neighbour == -1) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    vtkIdType id_face = m_Part.c2cGG(id_cell, i_face);
    if (id_face == -1) {
      EG_BUG;
    }
    if (!isSurface(id_cell, m_Grid)) {
      EG_BUG;
    }
    bc = cell_code->GetValue(id_face);
  }
  QList<node_t> nodes;
  nodes.append(node_t(id_node));
  nodes.append(node_t(id_node, edge_nodes[0]));
  nodes.append(node_t(face_nodes));
  nodes.append(node_t(id_node, edge_nodes[1]));
  vec3_t n = getNormalOfCell(m_Grid, id_cell, i_face);
  n.normalise();
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::createEdgeFace(vtkIdType id_node1, vtkIdType id_node2)
{
  // check if additional edge node needs to be created
  // (transition from boundary layer to far-field)
  // try to find a boundary layer cell which contains both nodes
  bool add_edge_node = false;
  for (int i = 0; i < m_Part.n2cGSize(id_node1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(id_node1, i);
    if (m_Cell2PCell[id_cell] != -1) {
      if (cellContainsNode(m_Grid, id_cell, id_node1) && cellContainsNode(m_Grid, id_cell, id_node2)) {
        add_edge_node = true;
        break;
      }
    }
  }
  if (!add_edge_node) {
    for (int i = 0; i < m_Part.n2cGSize(id_node2); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node2, i);
      if (m_Cell2PCell[id_cell] != -1) {
        if (cellContainsNode(m_Grid, id_cell, id_node1) && cellContainsNode(m_Grid, id_cell, id_node2)) {
          add_edge_node = true;
          break;
        }
      }
    }
  }

  QList<vtkIdType> cells;
  bool loop = false;
  getSortedEdgeCells(id_node1, id_node2, cells, loop);
  int owner     = m_Node2PCell[id_node1];
  int neighbour = m_Node2PCell[id_node2];
  if (owner == -1) {
    EG_BUG;
  }
  if (neighbour == -1) {
    EG_BUG;
  }
  if (owner > neighbour) {
    swap(id_node1, id_node2);
    swap(owner, neighbour);
  }
  QList<node_t> nodes;
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    node_t node;
    QVector<vtkIdType> pts1;
    getPointsOfCell(m_Grid, cells[i_cells], pts1);
    if (m_Cell2PCell[cells[i_cells]] != -1) {
      int i_cells2 = 0;
      if (i_cells == 0) {
        i_cells2 = 1;
      } else if (i_cells == cells.size() - 1) {
        i_cells2 = cells.size() - 2;
      } else {
        EG_BUG;
      }
      QVector<vtkIdType> pts2;
      getPointsOfCell(m_Grid, cells[i_cells2], pts2);
      QSet<vtkIdType> p1, p2;
      for (int i_pts1 = 0; i_pts1 < pts1.size(); ++i_pts1) {
        p1.insert(pts1[i_pts1]);
      }
      for (int i_pts2 = 0; i_pts2 < pts2.size(); ++i_pts2) {
        p2.insert(pts2[i_pts2]);
      }
      QSet<vtkIdType> face_nodes = p1.intersect(p2);
      node.id.resize(face_nodes.size());
      int i = 0;
      foreach (vtkIdType id_node, face_nodes) {
        node.id[i] = id_node;
        ++i;
      }
    } else {
      node.id.resize(pts1.size());
      for (int i_pts1 = 0; i_pts1 < pts1.size(); ++i_pts1) {
        node.id[i_pts1] = pts1[i_pts1];
      }
    }
    qSort(node.id);
    nodes.append(node);
  }
  if (!loop) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    if (cell_code->GetValue(cells.first()) != cell_code->GetValue(cells.last()) || add_edge_node) {
      nodes.append(node_t(id_node1, id_node2));
    }
  }
  vec3_t x1, x2;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  createFace(nodes, owner, neighbour, x2 - x1, 0);
}

void PolyMesh::createFaceFace(vtkIdType id_cell, int i_face)
{
  if (m_Cell2PCell[id_cell] == -1) {
    EG_BUG;
  }
  int owner = m_Cell2PCell[id_cell];
  int neighbour = -1;
  vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
  int bc = 0;
  if (id_neigh != -1) {
    if (isVolume(id_neigh, m_Grid)) {
      if (m_Cell2PCell[id_neigh] == -1) {
        EG_BUG;
      }
      neighbour = m_Cell2PCell[id_neigh];
    } else {
      EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
      bc = cell_code->GetValue(id_neigh);
    }
  }
  QVector<vtkIdType> tmp_node_ids;
  getFaceOfCell(m_Grid, id_cell, i_face, tmp_node_ids);
  vec3_t n = getNormalOfCell(m_Grid, id_cell, i_face);
  n.normalise();
  QVector<vtkIdType> node_ids(tmp_node_ids.size() + 1);
  for (int i = 0; i < tmp_node_ids.size(); ++i) {
    node_ids[i] = tmp_node_ids[i];
  }
  node_ids[node_ids.size() - 1] = node_ids[0];
  QList<node_t> nodes;
  for (int i = 0; i < node_ids.size() - 1; ++i) {
    nodes.append(node_t(node_ids[i]));
    if (m_Node2PCell[node_ids[i]] != -1 && m_Node2PCell[node_ids[i+1]] != -1) {
      nodes.append(node_t(node_ids[i], node_ids[i+1]));
    }
  }
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::createPointFace(vtkIdType id_node, int bc)
{
  bool is_loop;
  QList<vtkIdType> faces;
  getSortedPointFaces(id_node, bc, faces, is_loop);
  if (faces.size() == 0) {
    return;
  }
  QList<node_t> nodes;
  vec3_t n(0,0,0);
  if (faces.size() == 0) {
    EG_BUG;
  }
  foreach (vtkIdType id_face, faces) {
    node_t node;
    vtkIdType num_pts, *pts;
    m_Grid->GetCellPoints(id_face, num_pts, pts);
    node.id.resize(num_pts);
    for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
      node.id[i_pts] = pts[i_pts];
    }
    qSort(node.id);
    nodes.append(node);
    n += GeometryTools::cellNormal(m_Grid, id_face);
  }
  n.normalise();
  if (!is_loop) {
    bool prepend = false;
    vtkIdType id_face = faces.last();
    while (id_face != -1) {
      bool found = false;
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_face, num_pts, pts);
      QList<vtkIdType> id_neigh_node;
      vtkIdType id_neigh = -1;
      for (int i = 0; i < m_Part.c2cGSize(id_face); ++i) {
        id_neigh = m_Part.c2cGG(id_face, i);
        if (id_neigh == -1) {
          EG_BUG;
        }
        if (!isSurface(id_neigh, m_Grid)) {
          EG_BUG;
        }
        if (cellContainsNode(m_Grid, id_neigh, id_node)) {
          if (!faces.contains(id_neigh)) {
            if (found) {
              EG_BUG;
            }
            for (int j = 0; j < num_pts; ++j) {
              if (pts[j] != id_node && cellContainsNode(m_Grid, id_neigh, pts[j])) {
                id_neigh_node.append(pts[j]);
              }
            }
          }
        }
      }

      if (id_neigh_node.size() == 0) {
        EG_BUG;
      }

      if (id_neigh_node.size() == 1) {
        if (prepend) {
          nodes.prepend(node_t(id_node, id_neigh_node[0]));
          id_face = -1;
        } else {
          nodes.append(node_t(id_node, id_neigh_node[0]));
          prepend = true;
          id_face = faces.first();
        }
      } else {
        if (id_neigh_node.size() > 2) {
          EG_BUG;
        }
        if (faces.size() != 1) {
          EG_BUG;
        }
        nodes.prepend(node_t(id_node, id_neigh_node[0]));
        nodes.append(node_t(id_node, id_neigh_node[1]));
        id_face = -1;
      }
    }

    // check if this is a transition node (boundary layer -> far-field)
    bool dual_cell_found = false;
    bool cell_cell_found = false;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      if (m_Cell2PCell[m_Part.n2cGG(id_node,i)] == -1) {
        cell_cell_found = true;
      } else {
        dual_cell_found = true;
      }
    }

    if (m_Part.n2bcGSize(id_node) > 2 || (dual_cell_found && cell_cell_found)) {
      nodes.prepend(node_t(id_node));
    }
  }
  int owner     = m_Node2PCell[id_node];
  int neighbour = -1;
  createFace(nodes, owner, neighbour, n, bc);
}

void PolyMesh::computePoints()
{
  QVector<node_t> nodes;
  m_Nodes.getQVector(nodes);
  m_Points.resize(nodes.size());
  m_PointWeights.fill(1.0, m_Grid->GetNumberOfPoints());

  // find transition nodes
  QVector<bool> is_transition_node(m_Grid->GetNumberOfPoints(), false);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_Node2PCell[id_node] != -1) {
      for (int i_cell = 0; i_cell < m_Part.n2cGSize(id_node); ++i_cell) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i_cell);
        if (m_Cell2PCell[id_cell] != -1) {
          is_transition_node[id_node] = true;
          break;
        }
      }
    }
  }

  // compute weights (attraction for convex corners)
  // mark convex edge nodes
  QVector<bool> is_convex_node(m_Grid->GetNumberOfPoints(), false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Cell2PCell[id_cell] != -1) {
      vec3_t xc1 = cellCentre(m_Grid, id_cell);
      int N = 0;
      vec3_t n_out(0,0,0);
      for (int i_face = 0; i_face < m_Part.c2cGSize(id_cell); ++i_face) {
        vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
        if (id_neigh != -1) {
          if (m_Cell2PCell[id_neigh] == -1 && !isSurface(id_neigh, m_Grid)) {
            ++N;
            n_out += getNormalOfCell(m_Grid, id_cell, i_face);
          }
        }
      }
      if (N > 0) {
        n_out.normalise();
        for (int i_face = 0; i_face < m_Part.c2cGSize(id_cell); ++i_face) {
          vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
          if (id_neigh != -1) {
            if (m_Cell2PCell[id_neigh] != -1) {
              vec3_t xc2 = cellCentre(m_Grid, id_neigh);
              QVector<vtkIdType> face_nodes;
              getFaceOfCell(m_Grid, id_cell, i_face, face_nodes);
              if (face_nodes.size() == 0) {
                EG_BUG;
              }
              vec3_t xf(0,0,0);
              foreach (vtkIdType id_node, face_nodes) {
                vec3_t x;
                m_Grid->GetPoint(id_node, x.data());
                xf += x;
              }
              xf *= 1.0/face_nodes.size();
              vec3_t v1 = xc1 - xf;
              vec3_t v2 = xc2 - xf;
              v1.normalise();
              v2.normalise();
              if ((v1+v2)*n_out < 0) {
                double w = 1 + m_AttractorWeight*(1.0 + v1*v2);
                foreach (vtkIdType id_node, face_nodes) {
                  if (m_Node2PCell[id_node] != -1) {
                    m_PointWeights[id_node] = max(m_PointWeights[id_node], w);
                    if (v1*v2 > 0) {
                      //is_convex_node[id_node] = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // mapping for simple nodes
  QVector<int> prime2dual(m_Grid->GetNumberOfPoints(), -1);

  // compute the point locations
  for (int i = 0; i < nodes.size(); ++i) {
    if (nodes[i].id.size() == 0) {
      EG_BUG;
    }
    m_Points[i] = vec3_t(0,0,0);
    double weight_sum = 0.0;
    foreach (vtkIdType id, nodes[i].id) {
      vec3_t x;
      m_Grid->GetPoint(id, x.data());
      weight_sum += m_PointWeights[id];
      m_Points[i] += m_PointWeights[id]*x;
    }
    //m_Points[i] *= 1.0/nodes[i].id.size();
    m_Points[i] *= 1.0/weight_sum;
    if (nodes[i].id.size() == 1) {
      prime2dual[nodes[i].id[0]] = i;
    }
  }

  // correct transition nodes (pull-in)
  for (vtkIdType id_node1 = 0; id_node1 < m_Grid->GetNumberOfPoints(); ++id_node1) {
    if (is_transition_node[id_node1]) {
      if (prime2dual[id_node1] == -1) {
        EG_BUG;
      }
      vtkIdType id_node2 = -1;
      bool pull_in = false;
      for (int i_neigh = 0; i_neigh < m_Part.n2nGSize(id_node1); ++i_neigh) {
        vtkIdType id_neigh = m_Part.n2nGG(id_node1, i_neigh);
        if (m_Node2PCell[id_neigh] == -1 && !is_transition_node[id_neigh]) {
          if (id_node2 == -1) {
            pull_in = true;
            id_node2 = id_neigh;
          } else {
            pull_in = false;
          }
        }
      }
      if (pull_in) {
        double w = m_PullInFactor;
        m_Points[prime2dual[id_node1]] = w*m_Points[prime2dual[id_node2]] + (1-w)*m_Points[prime2dual[id_node1]];
      }
    }
  }

}

void PolyMesh::splitConcaveFaces()
{
  QList<int> delete_faces;
  do {
    delete_faces.clear();
    QList<face_t> new_faces;
    for (int i_face = 0; i_face < m_Faces.size(); ++i_face) {
      face_t face = m_Faces[i_face];
      int num_nodes = face.node.size();
      if (num_nodes >= 4) {
        QVector<vec3_t> x_face(num_nodes);
        for (int i = 0; i < num_nodes; ++i) {
          x_face[i] = nodeVector(face.node[i]);
        }
        vec3_t xc, n;
        GeometryTools::planeFit(x_face, xc, n);
        QVector<vec3_t> x(num_nodes + 2);
        for (int i = 0; i < num_nodes; ++i) {
          x[i+1] = x_face[i];
        }
        x.first() = x_face.last();
        x.last() = x_face.first();
        double L_max = 0.1;
        int i1 = -1;
        vec3_t v;
        for (int i = 1; i <= num_nodes; ++i) {
          vec3_t x_neigh = 0.5*(x[i-1] + x[i+1]);
          double Lnorm = (x[i-1] - x[i+1]).abs();
          double L = ((x_neigh - xc).abs() - (x[i] - xc).abs())/Lnorm;
          if (L > L_max) {
            i1 = i;
            L_max = L;
            v = x[i] - x_neigh;
          }
        }
        if (i1 != -1) {
          int i2 = -1;
          double alpha_min = 1e99;
          for (int i = 1; i <= num_nodes; ++i) {
            if ((i < i1 - 1 || i > i1 + 1) && (i < i1 + num_nodes -1 || i > i1 + num_nodes + 1)) {
              double alpha = GeometryTools::angle(x[i] - x[i1], v);
              if (alpha < alpha_min) {
                i2 = i;
                alpha_min = alpha;
              }
            }
          }
          if (i2 != -1) {
            --i1;
            --i2;
            if (i1 > i2) {
              swap(i1, i2);
            }
            delete_faces.append(i_face);
            QVector<int> nodes1;
            {
              int i = 0;
              while (i < num_nodes) {
                nodes1.append(face.node[i]);
                if (i == i1) {
                  i = i2;
                } else {
                  ++i;
                }
              }
            }
            QVector<int> nodes2;
            for (int i = i1; i <= i2; ++i) {
              nodes2.append(face.node[i]);
            }
            face_t face1(nodes1.size(), face.owner, face.neighbour, face.ref_vec, face.bc);
            face_t face2(nodes2.size(), face.owner, face.neighbour, face.ref_vec, face.bc);
            face1.node = nodes1;
            face2.node = nodes2;
            new_faces.append(face1);
            new_faces.append(face2);
          }
        }
      }
    }
    {
      int offset = 0;
      foreach (int i_face, delete_faces) {
        m_Faces.removeAt(i_face - offset);
        ++offset;
      }
    }
    m_Faces.append(new_faces);
    cout << delete_faces.size() << " concave faces have been split." << endl;
    delete_faces.clear();
  } while (delete_faces.size() > 0);
  sortFaces();
  buildPoint2Face();
  buildPCell2Face();
}

void PolyMesh::collectBoundaryConditions()
{
  QSet<int> bcs;
  foreach (face_t face, m_Faces) {
    if (face.bc != 0) {
      bcs.insert(face.bc);
    }
  }
  m_BCs.resize(bcs.size());
  qCopy(bcs.begin(), bcs.end(), m_BCs.begin());
}

void PolyMesh::createNodesAndFaces()
{
  m_Nodes.resize(m_Grid->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  int num_pcells = 0;
  int num_dcells = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Cell2PCell[id_cell] != -1) {
      ++num_pcells;
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_Node2PCell[id_node] != -1) {
      ++num_dcells;
    }
  }
  // create all remaining elements (prisms, hexes, and polyhedra)
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Cell2PCell[id_cell] != -1) {
      for (int i_face = 0; i_face < m_Part.c2cGSize(id_cell); ++i_face) {
        vtkIdType id_neigh = m_Part.c2cGG(id_cell, i_face);
        if (id_neigh == -1) {
          EG_BUG;
        }
        bool create_corner_faces = false;
        if (m_Cell2PCell[id_neigh] == -1 && !isSurface(id_neigh, m_Grid)) {
          create_corner_faces = true;
        }
        if (create_corner_faces) {
          QVector<vtkIdType> face_nodes;
          getFaceOfCell(m_Grid, id_cell, i_face, face_nodes);
          foreach (vtkIdType id_node, face_nodes) {
            if (m_Node2PCell[id_node] == -1) {
              EG_BUG;
            }
            createCornerFace(id_cell, i_face, id_node);
          }
        } else {
          if (id_neigh > id_cell || isSurface(id_neigh, m_Grid)) {
            createFaceFace(id_cell, i_face);
          }
        }
      }
    }
  }

  // create all dual cells
  if (m_CreateDualMesh) {
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      QSet<vtkIdType> c1;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (m_Cell2PCell[id_cell] == -1) {
          c1.insert(id_cell);
        }
      }
      if (m_Node2PCell[id_node] != -1) {
        for (int i_neigh = 0; i_neigh < m_Part.n2nGSize(id_node); ++i_neigh) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node, i_neigh);
          if (m_Node2PCell[id_neigh] != -1 && id_neigh > id_node) {

            // check if any of the adjacent cells (id_node <-> id_neigh) needs to be "dualised"
            QSet<vtkIdType> c2;
            for (int i = 0; i < m_Part.n2cGSize(id_neigh); ++i) {
              vtkIdType id_cell = m_Part.n2cGG(id_neigh, i);
              if (m_Cell2PCell[id_cell] == -1) {
                c2.insert(id_cell);
              }
            }
            c2.intersect(c1);
            if (c2.size() > 0) {
              createEdgeFace(id_node, id_neigh);
            }
          }
        }
        QSet<int> bcs;
        for (int i_cell = 0; i_cell < m_Part.n2cGSize(id_node); ++i_cell) {
          vtkIdType id_cell = m_Part.n2cGG(id_node, i_cell);
          if (isSurface(id_cell, m_Grid)) {
            bcs.insert(cell_code->GetValue(id_cell));
          }
        }
        foreach (int bc, bcs) {
          createPointFace(id_node, bc);
        }
      }
    }
  }

  computePoints();

  //splitConcaveFaces();
  sortFaces();
  collectBoundaryConditions();
}

vec3_t PolyMesh::faceNormal(int i)
{
  int N = m_Faces[i].node.size();
  QVector<vec3_t> x(N + 1);
  vec3_t xc(0,0,0);
  for (int j = 0; j < N; ++j) {
    x[j] = m_Points[m_Faces[i].node[j]];
    xc += x[j];
  }
  x[N] = x[0];
  xc *= 1.0/N;
  vec3_t n(0,0,0);
  for (int j = 0; j < N; ++j) {
    vec3_t u = x[j] - xc;
    vec3_t v = x[j+1] - xc;
    n += 0.5*u.cross(v);
  }
  return n;
}

void PolyMesh::checkFaceOrientation()
{
  for (int i = 0; i < m_Faces.size(); ++i) {
    int N = m_Faces[i].node.size();
    vec3_t n = faceNormal(i);
    n.normalise();
    if (n*m_Faces[i].ref_vec < 0) {
      QVector<int> old_node = m_Faces[i].node;
      for (int j = 0; j < N; ++j) {
        m_Faces[i].node[j] = old_node[N-1-j];
      }
    }
  }
}

void PolyMesh::buildPoint2Face()
{
  m_Point2Face.clear();
  m_Point2Face.resize(totalNumNodes());
  for (int face = 0; face < numFaces(); ++face) {
    foreach (int node, m_Faces[face].node) {
      m_Point2Face[node].append(face);
    }
  }
}

void PolyMesh::buildPCell2Face()
{
  m_PCell2Face.clear();
  m_PCell2Face.resize(m_NumPolyCells);
  for (int i = 0; i < m_Faces.size(); ++i) {
    if (owner(i) >= m_NumPolyCells) {
      EG_BUG;
    }
    m_PCell2Face[owner(i)].append(i);
    if (neighbour(i) >= 0) {
      if (neighbour(i) >= m_NumPolyCells) {
        EG_BUG;
      }
      m_PCell2Face[neighbour(i)].append(i);
    }
  }
}

void PolyMesh::invertFace(int i)
{
  invertQContainer(m_Faces[i].node);
  m_Faces[i].ref_vec *= -1;
}

void PolyMesh::sortFaces()
{
  // compute hashes
  int max_bc = -1;
  foreach (face_t face, m_Faces) {
    max_bc = max(max_bc, face.bc);
  }
  if (max_bc < 0) {
    EG_ERR_RETURN("mesh is corrupted");
  }
  int hash_stride = 1;
  foreach (face_t face, m_Faces) {
    hash_stride = max(hash_stride, face.owner);
  }
  QVector<QList<face_t> > sort_lists(hash_stride*(max_bc + 1) + 1);
  foreach (face_t face, m_Faces) {
    int i = face.bc*hash_stride + face.owner;
    sort_lists[i].append(face);
  }
  m_Faces.clear();
  for (int i = 0; i < sort_lists.size(); ++i) {
    qSort(sort_lists[i]);
    m_Faces += sort_lists[i];
  }
}
