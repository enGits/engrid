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

#include "createboundarylayershell.h"

#include "createvolumemesh.h"
#include "math/mathvector.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"
#include "deletevolumegrid.h"
#include "laplacesmoother.h"
#include "surfacemeshsmoother.h"
#include "deletecells.h"
#include "stitchholes.h"
#include "cgaltricadinterface.h"
#include "correctsurfaceorientation.h"
#include "geometrytools.h"
#include <vtkIdList.h>

namespace {

double triangleMinAngle(const vec3_t &a, const vec3_t &b, const vec3_t &c)
{
  vec3_t ab = b - a;
  vec3_t ac = c - a;
  vec3_t bc = c - b;
  if (ab.abs() < 1e-12 || ac.abs() < 1e-12 || bc.abs() < 1e-12) {
    return 0.0;
  }
  double ang_a = GeometryTools::angle(ab, ac);
  double ang_b = GeometryTools::angle(-1*ab, bc);
  double ang_c = GeometryTools::angle(-1*ac, -1*bc);
  return min(ang_a, min(ang_b, ang_c));
}

bool useQuadDiagonal02(const vec3_t &p0, const vec3_t &p1, const vec3_t &p2, const vec3_t &p3)
{
  double min_diag_02 = min(triangleMinAngle(p0, p1, p2), triangleMinAngle(p0, p2, p3));
  double min_diag_13 = min(triangleMinAngle(p0, p1, p3), triangleMinAngle(p1, p2, p3));
  return min_diag_02 >= min_diag_13;
}

void sanitizePolygonNodes(QVector<vtkIdType> &nodes)
{
  if (nodes.isEmpty()) {
    return;
  }
  QVector<vtkIdType> cleaned;
  cleaned.reserve(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    vtkIdType id = nodes[i];
    if (cleaned.isEmpty() || cleaned.last() != id) {
      cleaned.append(id);
    }
  }
  if (cleaned.size() > 1 && cleaned.first() == cleaned.last()) {
    cleaned.removeLast();
  }
  nodes = cleaned;
}

}


CreateBoundaryLayerShell::DeleteBadNodes::DeleteBadNodes(QList<vtkIdType> bad_nodes)
{
  m_BadNodes = bad_nodes;
  m_PerformGeometricChecks = false;
}

bool CreateBoundaryLayerShell::DeleteBadNodes::checkEdge(vtkIdType id_node1, vtkIdType id_node2)
{
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  if (m_BadNodes.contains(id_node1) && !m_BadNodes.contains(id_node2)) {
    return true;
  }
  return false;
}

CreateBoundaryLayerShell::CreateBoundaryLayerShell()
{
  m_MarchingDirection = 1;
  m_ReverseBoundaryFaces = false;
  m_RestGrid      = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_OriginalGrid  = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_PrismaticGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_Success = false;
}

void CreateBoundaryLayerShell::prepare()
{
  m_Part.trackGrid(m_Grid);

  DeleteVolumeGrid delete_volume;
  delete_volume.setGrid(m_Grid);
  delete_volume.setAllCells();
  delete_volume();

  m_Part.trackGrid(m_Grid);

  readSettings();
  setAllCells();
  getSurfaceCells(m_BoundaryLayerCodes, layer_cells, m_Grid);

  m_OriginalNodeNormals.resize(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_OriginalNodeNormals[id_node] = m_Part.globalNormal(id_node);
  }

  if (m_ReverseBoundaryFaces) {
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (m_BoundaryLayerCodes.contains(cell_code->GetValue(id_cell))) {
          reorientateFace(m_Grid, id_cell);
        }
      }
    }
  }

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (vtkIdType id_cell, layer_cells) {
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
      int bc = cell_code->GetValue(id_neigh);
      if (!m_BoundaryLayerCodes.contains(bc)) {
        m_LayerAdjacentBoundaryCodes.insert(bc);
      }
    }
  }

  // compute normals and origins of adjacent planes
  m_LayerAdjacentNormals.clear();
  m_LayerAdjacentOrigins.clear();
  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    double L = EG_LARGE_REAL;
    vec3_t n0(0, 0, 0);
    vec3_t x0(0, 0, 0);
    double total_area = 0;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid) && cell_code->GetValue(id_cell) == bc) {
        vec3_t n = cellNormal(m_Grid, id_cell);
        double A = n.abs();
        total_area += A;
        n0 += n;
        x0 += A*cellCentre(m_Grid, id_cell);
        L = min(L, sqrt(4*A/sqrt(3.0)));
      }
    }
    n0.normalise();
    x0 *= 1.0/total_area;
    m_LayerAdjacentNormals[bc] = n0;
    m_LayerAdjacentOrigins[bc] = x0;
  }

  computeBoundaryLayerVectors();
  makeCopy(m_Grid, m_OriginalGrid);
  if (m_PolyhedralBoundaryLayer) {
    dualiseBoundaryLayerSurface();
  }
  orientBoundaryLayerSurface();
}

void CreateBoundaryLayerShell::dualiseBoundaryLayerSurface()
{
  struct poly_t {
    QVector<vtkIdType> nodes;
    int bc;
  };

  if (!m_PolyhedralBoundaryLayer) {
    return;
  }

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  QVector<vec3_t> original_disp = m_BoundaryLayerVectors;
  QVector<vtkIdType> bl_faces;
  getSurfaceCells(m_BoundaryLayerCodes, bl_faces, m_Grid);
  if (bl_faces.isEmpty()) {
    return;
  }

  QSet<vtkIdType> bl_face_set;
  foreach (vtkIdType id_face, bl_faces) {
    bl_face_set.insert(id_face);
  }

  QMap<vtkIdType, QVector<vtkIdType> > face_nodes;
  QMap<vtkIdType, vec3_t> face_normal;
  QMap<vtkIdType, vtkIdType> face_centroid_node;

  for (int i = 0; i < bl_faces.size(); ++i) {
    vtkIdType id_face = bl_faces[i];
    EG_GET_CELL(id_face, m_Grid);
    QVector<vtkIdType> nodes(num_pts);
    for (int j = 0; j < num_pts; ++j) {
      nodes[j] = pts[j];
    }
    face_nodes[id_face] = nodes;
    vec3_t n = cellNormal(m_Grid, id_face);
    face_normal[id_face] = n;
  }

  QVector<vec3_t> dual_points;
  QVector<vec3_t> dual_vectors;
  m_DualNodeBoundaryCode.clear();
  m_DualNodeNormals.clear();

  auto add_dual_point = [&](const vec3_t &x, const vec3_t &d, const vec3_t &nref, int boundary_bc) -> vtkIdType {
    vtkIdType id = dual_points.size();
    dual_points.append(x);
    dual_vectors.append(d);
    m_DualNodeBoundaryCode.append(boundary_bc);
    vec3_t n = nref;
    if (n.abs() > 1e-12) {
      n.normalise();
    }
    m_DualNodeNormals.append(n);
    return id;
  };

  auto face_centroid = [&](vtkIdType id_face) -> vtkIdType {
    if (face_centroid_node.contains(id_face)) {
      return face_centroid_node[id_face];
    }
    EG_GET_CELL(id_face, m_Grid);
    vec3_t xc = cellCentre(m_Grid, id_face);
    vec3_t d(0,0,0);
    vec3_t nref(0,0,0);
    for (int i = 0; i < num_pts; ++i) {
      d += original_disp[pts[i]];
      if (pts[i] < m_OriginalNodeNormals.size()) {
        nref += m_OriginalNodeNormals[pts[i]];
      }
    }
    d *= 1.0/num_pts;
    vtkIdType id = add_dual_point(xc, d, nref, 0);
    face_centroid_node[id_face] = id;
    return id;
  };

  QMap<QPair<vtkIdType, vtkIdType>, vtkIdType> edge_midpoint_node;
  QMap<QPair<vtkIdType, vtkIdType>, int> edge_midpoint_bc;

  auto edge_key = [&](vtkIdType n1, vtkIdType n2) -> QPair<vtkIdType, vtkIdType> {
    if (n2 < n1) {
      qSwap(n1, n2);
    }
    return QPair<vtkIdType, vtkIdType>(n1, n2);
  };

  auto edge_midpoint = [&](vtkIdType n1, vtkIdType n2, int boundary_bc) -> vtkIdType {
    QPair<vtkIdType, vtkIdType> key = edge_key(n1, n2);
    if (edge_midpoint_node.contains(key)) {
      if (edge_midpoint_bc[key] == 0 && boundary_bc > 0) {
        edge_midpoint_bc[key] = boundary_bc;
        m_DualNodeBoundaryCode[edge_midpoint_node[key]] = boundary_bc;
      }
      return edge_midpoint_node[key];
    }
    vec3_t x1, x2;
    m_Grid->GetPoint(n1, x1.data());
    m_Grid->GetPoint(n2, x2.data());
    vec3_t x = 0.5*(x1 + x2);
    vec3_t d = 0.5*(original_disp[n1] + original_disp[n2]);
    vec3_t nref(0,0,0);
    if (n1 < m_OriginalNodeNormals.size()) {
      nref += m_OriginalNodeNormals[n1];
    }
    if (n2 < m_OriginalNodeNormals.size()) {
      nref += m_OriginalNodeNormals[n2];
    }
    vtkIdType id = add_dual_point(x, d, nref, boundary_bc);
    edge_midpoint_node[key] = id;
    edge_midpoint_bc[key] = boundary_bc;
    return id;
  };

  auto polygon_normal = [&](const QVector<vtkIdType> &poly_nodes) -> vec3_t {
    vec3_t xc(0,0,0);
    int N = poly_nodes.size();
    for (int i = 0; i < N; ++i) {
      xc += dual_points[poly_nodes[i]];
    }
    if (N > 0) {
      xc *= 1.0/N;
    }
    vec3_t n(0,0,0);
    for (int i = 0; i < N; ++i) {
      vec3_t u = dual_points[poly_nodes[i]] - xc;
      vec3_t v = dual_points[poly_nodes[(i + 1) % N]] - xc;
      n += 0.5*u.cross(v);
    }
    return n;
  };

  QVector<poly_t> polys;

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (!m_BoundaryLayerNode[id_node]) {
      continue;
    }
    QMap<int, QVector<vtkIdType> > faces_by_bc;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_face = m_Part.n2cGG(id_node, i);
      if (isSurface(id_face, m_Grid)) {
        int bc = cell_code->GetValue(id_face);
        if (m_BoundaryLayerCodes.contains(bc)) {
          faces_by_bc[bc].append(id_face);
        }
      }
    }

    QMap<int, QVector<vtkIdType> >::const_iterator bc_it;
    for (bc_it = faces_by_bc.constBegin(); bc_it != faces_by_bc.constEnd(); ++bc_it) {
      int bc = bc_it.key();
      QVector<vtkIdType> faces = bc_it.value();
      if (faces.isEmpty()) {
        continue;
      }

      QSet<vtkIdType> face_set;
      for (int i = 0; i < faces.size(); ++i) {
        face_set.insert(faces[i]);
      }

      QMap<vtkIdType, QSet<vtkIdType> > adj;
      QMap<vtkIdType, QVector<vtkIdType> > boundary_edges;
      QMap<QPair<vtkIdType, vtkIdType>, int> boundary_edge_bc;

      for (int i = 0; i < faces.size(); ++i) {
        vtkIdType id_face = faces[i];
        QVector<vtkIdType> nodes = face_nodes[id_face];
        int idx = nodes.indexOf(id_node);
        if (idx < 0) {
          continue;
        }
        vtkIdType prev = nodes[(idx + nodes.size() - 1) % nodes.size()];
        vtkIdType next = nodes[(idx + 1) % nodes.size()];
        vtkIdType edge_nodes[2] = { prev, next };
        for (int e = 0; e < 2; ++e) {
          vtkIdType id_other = edge_nodes[e];
          QSet<vtkIdType> edge_cells;
          getEdgeCells(id_node, id_other, edge_cells);
          QVector<vtkIdType> edge_bl_faces;
          foreach (vtkIdType id_edge_face, edge_cells) {
            if (isSurface(id_edge_face, m_Grid)) {
              if (cell_code->GetValue(id_edge_face) == bc && face_set.contains(id_edge_face)) {
                edge_bl_faces.append(id_edge_face);
              }
            }
          }

          bool boundary_edge = true;
          if (edge_bl_faces.size() == 2) {
            vtkIdType id_other_face = edge_bl_faces[0] == id_face ? edge_bl_faces[1] : edge_bl_faces[0];
            double ang = GeometryTools::angle(face_normal[id_face], face_normal[id_other_face]);
            if (ang < m_FeatureAngle) {
              adj[id_face].insert(id_other_face);
              adj[id_other_face].insert(id_face);
              boundary_edge = false;
            }
          }

          if (boundary_edge) {
            boundary_edges[id_face].append(id_other);
            int adj_bc = 0;
            foreach (vtkIdType id_edge_face, edge_cells) {
              if (isSurface(id_edge_face, m_Grid)) {
                int bc_edge = cell_code->GetValue(id_edge_face);
                if (!m_BoundaryLayerCodes.contains(bc_edge) && m_LayerAdjacentBoundaryCodes.contains(bc_edge)) {
                  adj_bc = bc_edge;
                  break;
                }
              }
            }
            QPair<vtkIdType, vtkIdType> key = edge_key(id_node, id_other);
            if (!boundary_edge_bc.contains(key) || boundary_edge_bc[key] == 0) {
              boundary_edge_bc[key] = adj_bc;
            }
          }
        }
      }

      QSet<vtkIdType> unvisited = face_set;
      while (!unvisited.isEmpty()) {
        vtkIdType start = *(unvisited.begin());
        QVector<vtkIdType> component;
        QList<vtkIdType> queue;
        queue.append(start);
        unvisited.remove(start);
        while (!queue.isEmpty()) {
          vtkIdType current = queue.takeFirst();
          component.append(current);
          foreach (vtkIdType neigh, adj[current]) {
            if (unvisited.contains(neigh)) {
              unvisited.remove(neigh);
              queue.append(neigh);
            }
          }
        }

        vtkIdType chain_start = -1;
        for (int i = 0; i < component.size(); ++i) {
          if (adj[component[i]].size() < 2) {
            chain_start = component[i];
            break;
          }
        }
        bool is_loop = chain_start == -1;
        if (is_loop) {
          chain_start = component[0];
        }

        QVector<vtkIdType> order;
        vtkIdType prev_face = -1;
        vtkIdType current_face = chain_start;
        while (true) {
          order.append(current_face);
          vtkIdType next_face = -1;
          foreach (vtkIdType neigh, adj[current_face]) {
            if (neigh != prev_face) {
              next_face = neigh;
              break;
            }
          }
          if (next_face == -1) {
            break;
          }
          prev_face = current_face;
          current_face = next_face;
          if (is_loop && current_face == chain_start) {
            break;
          }
        }

        QVector<vtkIdType> poly_nodes;

        auto boundary_other = [&](vtkIdType face_id, vtkIdType neigh_id) -> vtkIdType {
          QVector<vtkIdType> nodes = face_nodes[face_id];
          int idx = nodes.indexOf(id_node);
          if (idx < 0) {
            return -1;
          }
          vtkIdType prev = nodes[(idx + nodes.size() - 1) % nodes.size()];
          vtkIdType next = nodes[(idx + 1) % nodes.size()];
          if (neigh_id != -1) {
            QVector<vtkIdType> neigh_nodes = face_nodes[neigh_id];
            if (neigh_nodes.contains(prev)) {
              return next;
            }
            if (neigh_nodes.contains(next)) {
              return prev;
            }
          }
          if (boundary_edges.contains(face_id) && !boundary_edges[face_id].isEmpty()) {
            return boundary_edges[face_id].first();
          }
          return prev;
        };

        if (!is_loop) {
          if (order.size() == 1) {
            vtkIdType other1 = -1;
            vtkIdType other2 = -1;
            if (boundary_edges.contains(order[0]) && boundary_edges[order[0]].size() >= 2) {
              other1 = boundary_edges[order[0]][0];
              other2 = boundary_edges[order[0]][1];
          } else {
              other1 = boundary_other(order[0], -1);
              other2 = other1;
            }
            if (other1 < 0 || other2 < 0) {
              continue;
            }
            int bc1 = boundary_edge_bc.value(edge_key(id_node, other1), 0);
            int bc2 = boundary_edge_bc.value(edge_key(id_node, other2), 0);
            poly_nodes.append(edge_midpoint(id_node, other1, bc1));
            poly_nodes.append(face_centroid(order[0]));
            poly_nodes.append(edge_midpoint(id_node, other2, bc2));
          } else {
            vtkIdType start_face = order.first();
            vtkIdType end_face = order.last();
            vtkIdType start_other = boundary_other(start_face, order[1]);
            vtkIdType end_other = boundary_other(end_face, order[order.size() - 2]);
            if (start_other < 0 || end_other < 0) {
              continue;
            }
            int bc_start = boundary_edge_bc.value(edge_key(id_node, start_other), 0);
            int bc_end = boundary_edge_bc.value(edge_key(id_node, end_other), 0);
            poly_nodes.append(edge_midpoint(id_node, start_other, bc_start));
            for (int i = 0; i < order.size(); ++i) {
              poly_nodes.append(face_centroid(order[i]));
            }
            poly_nodes.append(edge_midpoint(id_node, end_other, bc_end));
          }
        } else {
          for (int i = 0; i < order.size(); ++i) {
            poly_nodes.append(face_centroid(order[i]));
          }
        }

        if (poly_nodes.size() < 3) {
          continue;
        }
        sanitizePolygonNodes(poly_nodes);
        if (poly_nodes.size() < 3) {
          continue;
        }

        vec3_t ref(0,0,0);
        for (int i = 0; i < poly_nodes.size(); ++i) {
          if (poly_nodes[i] < m_DualNodeNormals.size()) {
            ref += m_DualNodeNormals[poly_nodes[i]];
          }
        }
        if (ref.abs() <= 1e-12 && id_node < m_OriginalNodeNormals.size()) {
          ref = m_OriginalNodeNormals[id_node];
        }
        if (ref.abs() <= 1e-12) {
          for (int i = 0; i < order.size(); ++i) {
            ref += face_normal[order[i]];
          }
        }
        vec3_t n_poly = polygon_normal(poly_nodes);
        if (ref.abs() > 1e-12 && n_poly.abs() > 1e-12 && n_poly*ref < 0) {
          invertQContainer(poly_nodes);
        }

        poly_t poly;
        poly.nodes = poly_nodes;
        poly.bc = bc;
        polys.append(poly);
      }
    }
  }

  EG_VTKSP(vtkUnstructuredGrid, dual_grid);
  if (polys.isEmpty()) {
    m_DualNodeBoundaryCode.clear();
    return;
  }
  allocateGrid(dual_grid, polys.size(), dual_points.size());
  for (vtkIdType id_node = 0; id_node < dual_points.size(); ++id_node) {
    dual_grid->GetPoints()->SetPoint(id_node, dual_points[id_node].data());
  }
  EG_VTKDCC(vtkIntArray, dual_cell_code, dual_grid, "cell_code");
  for (int i = 0; i < polys.size(); ++i) {
    const poly_t &poly = polys[i];
    EG_VTKSP(vtkIdList, pts);
    pts->SetNumberOfIds(poly.nodes.size());
    for (int j = 0; j < poly.nodes.size(); ++j) {
      pts->SetId(j, poly.nodes[j]);
    }
    vtkIdType id_new_cell = dual_grid->InsertNextCell(VTK_POLYGON, pts);
    dual_cell_code->SetValue(id_new_cell, poly.bc);
  }

  makeCopy(dual_grid, m_Grid);
  m_Grid->Modified();
  updateCellIndex(m_Grid);

  if (m_Grid->GetNumberOfCells() > 0) {
    CorrectSurfaceOrientation surf_check;
    surf_check.setGrid(m_Grid);
    surf_check.setStart(0);
    surf_check();
  }

  m_BoundaryLayerVectors = dual_vectors;
  m_Height.resize(dual_vectors.size());
  m_BoundaryLayerNode.fill(true, m_Grid->GetNumberOfPoints());
  for (int i = 0; i < dual_vectors.size(); ++i) {
    m_Height[i] = dual_vectors[i].abs();
    if (m_Height[i] <= 1e-12) {
      m_BoundaryLayerVectors[i] = vec3_t(0,0,0);
    }
  }

  m_Part.trackGrid(m_Grid);
}

void CreateBoundaryLayerShell::orientBoundaryLayerSurface()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, m_Grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  QVector<vtkIdType> bl_faces;
  getSurfaceCells(m_BoundaryLayerCodes, bl_faces, m_Grid);
  if (bl_faces.isEmpty()) {
    return;
  }

  foreach (vtkIdType id_face, bl_faces) {
    EG_GET_CELL(id_face, m_Grid);
    bool is_bl = m_BoundaryLayerCodes.contains(cell_code->GetValue(id_face));
    vec3_t ref(0,0,0);
    for (int i = 0; i < num_pts; ++i) {
      vtkIdType id_node = pts[i];
      if (m_PolyhedralBoundaryLayer) {
        if (id_node < m_DualNodeNormals.size()) {
          ref += m_DualNodeNormals[id_node];
        }
      } else {
        if (id_node < m_OriginalNodeNormals.size()) {
          ref += m_OriginalNodeNormals[id_node];
        }
      }
    }
    if (ref.abs() <= 1e-12) {
      continue;
    }
    if (m_ReverseBoundaryFaces && is_bl) {
      ref *= -1.0;
    }
    vec3_t n = cellNormal(m_Grid, id_face);
    if (n.abs() <= 1e-12) {
      continue;
    }
    if (n*ref < 0) {
      reorientateFace(m_Grid, id_face);
    }
    cell_orgdir->SetValue(id_face, cell_curdir->GetValue(id_face));
  }
}

QList<vtkIdType> CreateBoundaryLayerShell::findBadNodes(int bc)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  CgalTriCadInterface cad(m_ShellGrid);
  QList<vtkIdType> bad_nodes;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (!m_BoundaryLayerNode[id_node]){
      bool node_has_bc = false;
      bool all_bcs_adjacent = true;
      for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
        int bcn = m_Part.n2bcG(id_node, i);
        if (bcn == bc) {
          node_has_bc = true;
        }
        if (!m_LayerAdjacentBoundaryCodes.contains(bcn)) {
          all_bcs_adjacent = false;
        }
      }
      if (node_has_bc && all_bcs_adjacent && node_type->GetValue(id_node) != EG_FIXED_VERTEX) {
        vec3_t x1;
        m_Grid->GetPoint(id_node, x1.data());
        vec3_t x2 = cad.snap(x1);
        if ((x2-x1)*cad.getLastNormal() <= 0) {
          bad_nodes << id_node;
        }
      }
    }
  }
  return bad_nodes;
}

QList<vtkIdType> CreateBoundaryLayerShell::correctAdjacentBC(int bc, int num_levels)
{
  cout << "correcting boundary \"" << qPrintable(GuiMainWindow::pointer()->getBC(bc).getName()) << "\" with " << num_levels << " levels" << endl;

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  double scal_min = -1;
  int    count    = 0;

  SurfaceMeshSmoother smooth;
  smooth.setGrid(m_Grid);
  smooth.useSimpleCentreScheme();

  // mark num_levels levels of nodes next to the boundary layer
  QVector<bool> marked(m_Grid->GetNumberOfPoints(), false);
  QList<vtkIdType> marked_nodes;
  EG_FORALL_NODES(id_node, m_Grid) {
    if (m_BoundaryLayerNode[id_node]) {
      marked[id_node] = true;
      marked_nodes << id_node;
    }
  }
  for (int level = 0; level < num_levels; ++level) {
    QList<vtkIdType> new_marked_nodes;
    foreach (vtkIdType id_node, marked_nodes) {
      for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
        vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
        if (!marked[id_neigh]) {
          marked[id_neigh] = true;
          new_marked_nodes << id_neigh;
        }
      }
    }
    marked_nodes = new_marked_nodes;
  }

  QList<vtkIdType> cells;
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid)) {
      if (cell_code->GetValue(id_cell) == bc) {
        cells << id_cell;
      }
    }
  }

  CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(bc);
  smooth.setCells(cells);
  smooth.prepareCadInterface(cad);

  QList<vtkIdType> bad_nodes;
  int new_num_bad = m_Grid->GetNumberOfPoints();
  int old_num_bad = 0;
  updateNodeInfo();
  do {
    old_num_bad = new_num_bad;
    scal_min = 1;
    bad_nodes.clear();
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!m_BoundaryLayerNode[id_node]){
        bool node_has_bc = false;
        bool all_bcs_adjacent = true;
        for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
          int bcn = m_Part.n2bcG(id_node, i);
          if (bcn == bc) {
            node_has_bc = true;
          }
          if (!m_LayerAdjacentBoundaryCodes.contains(bcn)) {
            all_bcs_adjacent = false;
          }
        }
        if (node_has_bc && all_bcs_adjacent && marked[id_node] && node_type->GetValue(id_node) != EG_FIXED_VERTEX) {

          // one boundary condition and EG_SIMPLE_VERTEX
          //
          if (m_Part.n2bcGSize(id_node) == 1) {
            vec3_t xs(0, 0, 0); // = smooth.smoothNode(id_node);
            int count = 0;
            for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
              vtkIdType id_cell = m_Part.n2cGG(id_node, i);
              if (cell_code->GetValue(id_cell) == bc) {
                vec3_t x = cellCentre(m_Grid, id_cell);
                xs += x;
                ++count;
              }
            }
            if (count > 0) {
              xs *= 1.0/count;
              if (node_type->GetValue(id_node) == EG_SIMPLE_VERTEX) {
                vec3_t n = m_Part.globalNormal(id_node);
                if (m_OriginalNodeNormals[id_node]*n < 0) {
                  n = m_OriginalNodeNormals[id_node];
                }
                xs = cad->snapWithNormal(xs, m_OriginalNodeNormals[id_node]);
              } else {
                xs = cad->snapToEdge(xs);
              }
              if (!checkVector(xs)) {
                EG_ERR_RETURN("error while correcting adjacent boundaries");
              }
              m_Grid->GetPoints()->SetPoint(id_node, xs.data());
            }
          }

          // one boundary condition and EG_FEATURE_EDGE_VERTEX
          //
          else if (m_Part.n2bcGSize(id_node) == 1 && node_type->GetValue(id_node) == EG_FEATURE_EDGE_VERTEX) {
            int count = 0;
            vec3_t xs(0, 0, 0);
            for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
              vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
              if (node_type->GetValue(id_neigh) == EG_FEATURE_EDGE_VERTEX) {
                vec3_t x;
                m_Grid->GetPoint(id_neigh, x.data());
                xs += x;
                ++count;
              }
            }
            if (count < 2) {
              vec3_t x;
              m_Grid->GetPoint(id_node, x.data());
              xs += x;
              ++count;
            }
            xs *= 1.0/count;
            xs = cad->snapToEdge(xs);
            m_Grid->GetPoints()->SetPoint(id_node, xs.data());
          }

          // two boundary conditions and no EG_FIXED_VERTEX
          //
          else if (m_Part.n2bcGSize(id_node) == 2 && node_type->GetValue(id_node) != EG_FIXED_VERTEX) {
            int count = 0;
            vec3_t xs(0, 0, 0);
            for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
              vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
              if (m_Part.n2bcGSize(id_neigh) > 1) {
                vec3_t x;
                m_Grid->GetPoint(id_neigh, x.data());
                xs += x;
                ++count;
              }
            }
            if (count < 2) {
              vec3_t x;
              m_Grid->GetPoint(id_node, x.data());
              xs += x;
              ++count;
            }
            xs *= 1.0/count;
            xs = cad->snapToEdge(xs);
            m_Grid->GetPoints()->SetPoint(id_node, xs.data());
          }

          bool node_bad = false;
          for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
            vtkIdType id_cell = m_Part.n2cGG(id_node, i);
            if (isSurface(id_cell, m_Grid)) {
              if (cell_code->GetValue(id_cell) == bc) {
                CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(cell_code->GetValue(id_cell));
                cad->snap(cellCentre(m_Grid, id_cell));
                vec3_t n = cellNormal(m_Grid, id_cell);
                n.normalise();
                double scal = n*cad->getLastNormal();
                scal_min = min(scal_min, scal);
                if (scal < 0.5 && !node_bad) {
                  bad_nodes << id_node;
                  node_bad = true;
                }
              }
            }
          }
        }
      }
    }
    new_num_bad = bad_nodes.size();
    ++count;
  } while (scal_min < 0.5 && count < 20);
  cout << "  " << bad_nodes.size() << " node defects" << endl;
  return bad_nodes;
}

void CreateBoundaryLayerShell::createLayerNodes(vtkIdType id_node)
{
  vec3_t x1;
  m_Grid->GetPoint(id_node, x1.data());
  m_PrismaticGrid->GetPoints()->SetPoint(m_ShellNodeMap[id_node], x1.data());
  vec3_t x2 = x1 + m_MarchingDirection*m_BoundaryLayerVectors[id_node];

  m_PrismaticGrid->GetPoints()->SetPoint(m_ShellNodeMap[id_node], x1.data());
  double h1 = 0;
  for (int i = 1; i < m_NumLayers; ++i) {
    vec3_t x = x1 + m_RelativeHeights[i]*(x2 - x1);
    m_PrismaticGrid->GetPoints()->SetPoint(i*m_ShellPart.getNumberOfNodes() + m_ShellNodeMap[id_node], x.data());
  }
  m_PrismaticGrid->GetPoints()->SetPoint(m_NumLayers*m_ShellPart.getNumberOfNodes() + m_ShellNodeMap[id_node], x2.data());
  m_Grid->GetPoints()->SetPoint(id_node, x2.data());
}

void CreateBoundaryLayerShell::createPrismaticGrid()
{
  QVector<vtkIdType> original_faces, shell_faces;
  getSurfaceCells(m_BoundaryLayerCodes, original_faces, m_OriginalGrid);
  getSurfaceCells(m_BoundaryLayerCodes, shell_faces, m_Grid);
  vtkIdType extra_nodes = 0;
  if (m_PolyhedralBoundaryLayer && m_TriangulateFarfield) {
    foreach (vtkIdType id_cell, shell_faces) {
      EG_GET_CELL(id_cell, m_Grid);
      if (num_pts >= 5) {
        ++extra_nodes;
      }
    }
  }
  {
    MeshPartition part(m_Grid);
    part.setCells(shell_faces);
    allocateGrid(m_PrismaticGrid,
                 (m_NumLayers + 1)*part.getNumberOfCells(),
                 (m_NumLayers + 1)*part.getNumberOfNodes() + extra_nodes);
  }

  EG_VTKDCC(vtkIntArray, cell_code_pri, m_PrismaticGrid, "cell_code");

  m_ShellNodeMap.fill(-1, m_Grid->GetNumberOfPoints());
  m_ShellPart.setGrid(m_Grid);
  m_ShellPart.setCells(shell_faces);
  vtkIdType next_extra_node = (m_NumLayers + 1)*m_ShellPart.getNumberOfNodes();
  for (int i = 0; i < m_ShellPart.getNumberOfNodes(); ++i) {
    m_ShellNodeMap[m_ShellPart.globalNode(i)] = i;
  }

  QVector<QSet<int> > n2bc(m_PrismaticGrid->GetNumberOfPoints());

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_ShellNodeMap[id_node] != -1) {
      createLayerNodes(id_node);
      if (!m_PolyhedralBoundaryLayer) {
        for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
          n2bc[m_ShellNodeMap[id_node]].insert(m_Part.n2bcG(id_node, i));
        }
      }
    }
  }

  QList<QVector<vtkIdType> > adjacent_edges;

  // create prismatic cells and prepare adjacent quad faces
  //
  auto edge_boundary_code = [&](vtkIdType n1, vtkIdType n2) -> int {
    int bc1 = 0;
    int bc2 = 0;
    if (n1 >= 0 && n1 < m_DualNodeBoundaryCode.size() && m_DualNodeBoundaryCode[n1] > 0) {
      bc1 = m_DualNodeBoundaryCode[n1];
    }
    if (n2 >= 0 && n2 < m_DualNodeBoundaryCode.size() && m_DualNodeBoundaryCode[n2] > 0) {
      bc2 = m_DualNodeBoundaryCode[n2];
    }
    if (bc1 > 0 && bc2 > 0) {
      return bc1 == bc2 ? bc1 : bc1;
    }
    return bc1 > 0 ? bc1 : bc2;
  };

  foreach (vtkIdType id_cell, shell_faces) {
    EG_GET_CELL(id_cell, m_Grid);
    QVector<vtkIdType> face_pts(num_pts);
    for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
      if (m_ShellNodeMap[pts[i_pts]] < 0) {
        EG_BUG;
      }
      if (m_ShellNodeMap[pts[i_pts]] >= m_ShellPart.getNumberOfNodes()) {
        EG_BUG;
      }
      face_pts[i_pts] = m_ShellNodeMap[pts[i_pts]];
    }

    for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
      int i_next = (i_pts + 1) % num_pts;
      QVector<vtkIdType> edge(4);
      edge[1] = face_pts[i_pts];
      edge[2] = face_pts[i_next];
      int bc_edge = 0;
      if (m_PolyhedralBoundaryLayer) {
        bc_edge = edge_boundary_code(edge[1], edge[2]);
      } else {
        QSet<int> edge_codes = m_LayerAdjacentBoundaryCodes;
        edge_codes.intersect(n2bc[edge[1]]);
        edge_codes.intersect(n2bc[edge[2]]);
        if (edge_codes.size() == 1) {
          bc_edge = *edge_codes.begin();
        }
      }
      if (bc_edge > 0) {
        edge[0] = bc_edge;
        edge[3] = id_cell;
        adjacent_edges.append(edge);
      }
    }

    if (num_pts == 3) {
      vtkIdType tri_pts[3];
      for (int i = 0; i < 3; ++i) {
        tri_pts[i] = face_pts[i];
      }
      vtkIdType id_tri = m_PrismaticGrid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
      copyCellData(m_Grid, id_cell, m_PrismaticGrid, id_tri);
    } else if (num_pts == 4) {
      vtkIdType quad_pts[4];
      for (int i = 0; i < 4; ++i) {
        quad_pts[i] = face_pts[i];
      }
      vtkIdType id_quad = m_PrismaticGrid->InsertNextCell(VTK_QUAD, 4, quad_pts);
      copyCellData(m_Grid, id_cell, m_PrismaticGrid, id_quad);
    } else {
      EG_VTKSP(vtkIdList, poly_pts);
      poly_pts->SetNumberOfIds(num_pts);
      for (int i = 0; i < num_pts; ++i) {
        poly_pts->SetId(i, face_pts[i]);
      }
      vtkIdType id_poly = m_PrismaticGrid->InsertNextCell(VTK_POLYGON, poly_pts);
      copyCellData(m_Grid, id_cell, m_PrismaticGrid, id_poly);
    }

    for (int i_layer = 0; i_layer < m_NumLayers; ++i_layer) {
      bool triangulate_top = m_PolyhedralBoundaryLayer && m_TriangulateFarfield
                             && i_layer == m_NumLayers - 1 && num_pts > 3;
      bool make_polyhedron = num_pts > 4;
      if (triangulate_top) {
        make_polyhedron = true;
      }
      if (!make_polyhedron && num_pts == 3) {
        vtkIdType pri_pts[6];
        for (int i_pts = 0; i_pts < 3; ++i_pts) {
          pri_pts[i_pts]     = face_pts[i_pts] + i_layer*m_ShellPart.getNumberOfNodes();
          pri_pts[i_pts + 3] = face_pts[i_pts] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
        }
        vtkIdType id_pri = m_PrismaticGrid->InsertNextCell(VTK_WEDGE, 6, pri_pts);
        cell_code_pri->SetValue(id_pri, 0);
      } else if (!make_polyhedron && num_pts == 4) {
        vtkIdType pri_pts[8];
        for (int i_pts = 0; i_pts < 4; ++i_pts) {
          pri_pts[i_pts]     = face_pts[i_pts] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
          pri_pts[i_pts + 4] = face_pts[i_pts] + i_layer*m_ShellPart.getNumberOfNodes();
        }
        vtkIdType id_hex = m_PrismaticGrid->InsertNextCell(VTK_HEXAHEDRON, 8, pri_pts);
        cell_code_pri->SetValue(id_hex, 0);
      } else {
        vtkIdType bottom_offset = i_layer*m_ShellPart.getNumberOfNodes();
        vtkIdType top_offset = (i_layer + 1)*m_ShellPart.getNumberOfNodes();
        int top_faces = 1;
        if (triangulate_top) {
          if (num_pts == 4) {
            top_faces = 2;
          } else {
            top_faces = num_pts;
          }
        }
        int num_faces = 1 + num_pts + top_faces;
        QVector<vtkIdType> stream_ids;
        stream_ids.reserve(1 + num_faces*(num_pts + 2));
        stream_ids.append(num_faces);

        // bottom face (outward)
        stream_ids.append(num_pts);
        for (int i = 0; i < num_pts; ++i) {
          stream_ids.append(face_pts[i] + bottom_offset);
        }

        // top face(s)
        if (!triangulate_top) {
          stream_ids.append(num_pts);
          for (int i = num_pts - 1; i >= 0; --i) {
            stream_ids.append(face_pts[i] + top_offset);
          }
        } else if (num_pts == 4) {
          vec3_t p[4];
          for (int i = 0; i < 4; ++i) {
            m_PrismaticGrid->GetPoint(face_pts[i] + top_offset, p[i].data());
          }
          bool use_diag_02 = useQuadDiagonal02(p[0], p[1], p[2], p[3]);
          if (use_diag_02) {
            stream_ids.append(3);
            stream_ids.append(face_pts[0] + top_offset);
            stream_ids.append(face_pts[2] + top_offset);
            stream_ids.append(face_pts[1] + top_offset);
            stream_ids.append(3);
            stream_ids.append(face_pts[0] + top_offset);
            stream_ids.append(face_pts[3] + top_offset);
            stream_ids.append(face_pts[2] + top_offset);
          } else {
            stream_ids.append(3);
            stream_ids.append(face_pts[0] + top_offset);
            stream_ids.append(face_pts[3] + top_offset);
            stream_ids.append(face_pts[1] + top_offset);
            stream_ids.append(3);
            stream_ids.append(face_pts[1] + top_offset);
            stream_ids.append(face_pts[3] + top_offset);
            stream_ids.append(face_pts[2] + top_offset);
          }
        } else {
          vtkIdType center_id = next_extra_node++;
          vec3_t center(0, 0, 0);
          for (int i = 0; i < num_pts; ++i) {
            vec3_t x;
            m_PrismaticGrid->GetPoint(face_pts[i] + top_offset, x.data());
            center += x;
          }
          center *= (1.0/num_pts);
          m_PrismaticGrid->GetPoints()->SetPoint(center_id, center.data());
          for (int i = 0; i < num_pts; ++i) {
            vtkIdType p1 = face_pts[i] + top_offset;
            vtkIdType p2 = face_pts[(i + 1) % num_pts] + top_offset;
            stream_ids.append(3);
            stream_ids.append(p1);
            stream_ids.append(center_id);
            stream_ids.append(p2);
          }
        }

        // side faces
        for (int i = 0; i < num_pts; ++i) {
          stream_ids.append(4);
          vtkIdType p1 = face_pts[i] + bottom_offset;
          vtkIdType p2 = face_pts[(i + 1) % num_pts] + bottom_offset;
          vtkIdType p3 = face_pts[(i + 1) % num_pts] + top_offset;
          vtkIdType p4 = face_pts[i] + top_offset;
          stream_ids.append(p1);
          stream_ids.append(p4);
          stream_ids.append(p3);
          stream_ids.append(p2);
        }

        EG_VTKSP(vtkIdList, stream);
        stream->SetNumberOfIds(stream_ids.size());
        for (int i = 0; i < stream_ids.size(); ++i) {
          stream->SetId(i, stream_ids[i]);
        }
        vtkIdType id_poly = m_PrismaticGrid->InsertNextCell(VTK_POLYHEDRON, stream);
        cell_code_pri->SetValue(id_poly, 0);
      }
    }
  }

  // create quads on adjacent boundary faces
  //
  EG_VTKSP(vtkUnstructuredGrid, noquad_grid);
  makeCopy(m_PrismaticGrid, noquad_grid);
  allocateGrid(m_PrismaticGrid, m_PrismaticGrid->GetNumberOfCells() + m_NumLayers*adjacent_edges.size(), m_PrismaticGrid->GetNumberOfPoints());
  makeCopyNoAlloc(noquad_grid, m_PrismaticGrid);

  EG_VTKDCC(vtkIntArray, cell_code, m_PrismaticGrid, "cell_code");

  foreach (QVector<vtkIdType> edge, adjacent_edges) {
    vtkIdType qua_pts[4];
    for (int i_layer = 0; i_layer < m_NumLayers; ++i_layer) {
      qua_pts[0] = edge[2] + i_layer*m_ShellPart.getNumberOfNodes();
      qua_pts[1] = edge[1] + i_layer*m_ShellPart.getNumberOfNodes();
      qua_pts[2] = edge[1] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
      qua_pts[3] = edge[2] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
      vtkIdType id_qua = m_PrismaticGrid->InsertNextCell(VTK_QUAD, 4, qua_pts);
      copyCellData(m_Grid, edge[3], m_PrismaticGrid, id_qua);
      cell_code->SetValue(id_qua, edge[0]);
    }
  }
}

void CreateBoundaryLayerShell::triangulateFarfieldSurface(vtkUnstructuredGrid *shell_grid)
{
  if (!shell_grid || shell_grid->GetNumberOfCells() == 0) {
    return;
  }

  vtkIdType num_tri_cells = 0;
  vtkIdType extra_nodes = 0;
  for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
    if (!isSurface(id_cell, shell_grid)) {
      continue;
    }
    EG_GET_CELL(id_cell, shell_grid);
    if (num_pts < 3) {
      continue;
    }
    if (num_pts == 3) {
      ++num_tri_cells;
    } else if (num_pts == 4) {
      num_tri_cells += 2;
    } else {
      num_tri_cells += num_pts;
      ++extra_nodes;
    }
  }

  if (num_tri_cells == 0) {
    return;
  }

  EG_VTKSP(vtkUnstructuredGrid, tri_grid);
  allocateGrid(tri_grid, num_tri_cells, shell_grid->GetNumberOfPoints() + extra_nodes);
  for (vtkIdType id_node = 0; id_node < shell_grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    shell_grid->GetPoint(id_node, x.data());
    tri_grid->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(shell_grid, id_node, tri_grid, id_node);
  }

  vtkIdType next_extra_node = shell_grid->GetNumberOfPoints();
  for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
    if (!isSurface(id_cell, shell_grid)) {
      continue;
    }
    EG_GET_CELL(id_cell, shell_grid);
    if (num_pts < 3) {
      continue;
    }
    if (num_pts == 3) {
      vtkIdType tri_pts[3] = { pts[0], pts[1], pts[2] };
      vtkIdType id_new = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
      copyCellData(shell_grid, id_cell, tri_grid, id_new);
    } else if (num_pts == 4) {
      vec3_t p[4];
      for (int i = 0; i < 4; ++i) {
        shell_grid->GetPoint(pts[i], p[i].data());
      }
      bool use_diag_02 = useQuadDiagonal02(p[0], p[1], p[2], p[3]);
      if (use_diag_02) {
        vtkIdType tri_pts1[3] = { pts[0], pts[1], pts[2] };
        vtkIdType tri_pts2[3] = { pts[0], pts[2], pts[3] };
        vtkIdType id_new1 = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts1);
        vtkIdType id_new2 = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts2);
        copyCellData(shell_grid, id_cell, tri_grid, id_new1);
        copyCellData(shell_grid, id_cell, tri_grid, id_new2);
      } else {
        vtkIdType tri_pts1[3] = { pts[0], pts[1], pts[3] };
        vtkIdType tri_pts2[3] = { pts[1], pts[2], pts[3] };
        vtkIdType id_new1 = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts1);
        vtkIdType id_new2 = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts2);
        copyCellData(shell_grid, id_cell, tri_grid, id_new1);
        copyCellData(shell_grid, id_cell, tri_grid, id_new2);
      }
    } else {
      vtkIdType center_id = next_extra_node++;
      vec3_t center(0, 0, 0);
      for (int i = 0; i < num_pts; ++i) {
        vec3_t x;
        shell_grid->GetPoint(pts[i], x.data());
        center += x;
      }
      center *= (1.0/num_pts);
      tri_grid->GetPoints()->SetPoint(center_id, center.data());
      for (int i = 0; i < num_pts; ++i) {
        vtkIdType tri_pts[3] = { pts[i], pts[(i + 1) % num_pts], center_id };
        vtkIdType id_new = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
        copyCellData(shell_grid, id_cell, tri_grid, id_new);
      }
    }
  }

  makeCopy(tri_grid, shell_grid);
}

void CreateBoundaryLayerShell::reduceSurface()
{
  RemovePoints remove_points;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  remove_points.setMeshPartition(part);
  remove_points.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  remove_points.setUpdatePSPOn();
  remove_points.setThreshold(3);
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), false);
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    // reset all node types to EG_SIMPLE_VERTEX to trigger a recalculation of the node types
    node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_QUAD) {
        fix[id_node] = true;
      }
    }
  }
  remove_points.fixNodes(fix);
  remove_points();
}

void CreateBoundaryLayerShell::smoothSurface()
{
  LaplaceSmoother smooth;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  smooth.setMeshPartition(part);
  smooth.setNumberOfIterations(2);
  smooth.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), false);
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    // reset all node types to EG_SIMPLE_VERTEX to trigger a recalculation of the node types
    node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_QUAD) {
        fix[id_node] = true;
      }
    }
  }
  smooth.fixNodes(fix);
  smooth();
}

void CreateBoundaryLayerShell::operate()
{
  prepare();
  createPrismaticGrid();
  m_Success = true;
  m_Part.trackGrid(m_Grid);
  return;

  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    QList<vtkIdType> bad_nodes;
    int levels = 1;
    QVector<vec3_t> x_save(m_Grid->GetNumberOfPoints());
    do {
      ++levels;
      for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
        m_Grid->GetPoint(id_node, x_save[id_node].data());
      }
      QList<vtkIdType> last_bad = bad_nodes;
      bad_nodes = correctAdjacentBC(bc, levels);
      if (bad_nodes.size() >= last_bad.size() && last_bad.size() > 0) {
        bad_nodes = last_bad;
        for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
          m_Grid->GetPoints()->SetPoint(id_node, x_save[id_node].data());
        }
        break;
      }
    } while (bad_nodes.size() > 0);
  }

  QList<vtkIdType> bad_cells;
  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    QList<vtkIdType> bad_nodes = findBadNodes(bc);
    if (bad_nodes.size() > 0) {
      bool fixable = true;
      QVector<bool> is_bad_cell(m_Grid->GetNumberOfCells(), false);
      foreach (vtkIdType id_node, bad_nodes) {
        if (m_Part.n2bcGSize(id_node) != 1) {
          cout << "node " << id_node << " cannot be fixed." << endl;
          fixable = false;
          break;
        }
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          is_bad_cell[m_Part.n2cGG(id_node, i)] = true;
        }
      }
      if (fixable) {
        EG_FORALL_CELLS (id_cell, m_Grid) {
          if (is_bad_cell[id_cell]) {
            bad_cells << id_cell;
          }
        }
      } else {
        m_Success = false;
        cout << "adjacent patch cannot be corrected!" << endl;
        return;
      }
    }
  }

  cout << "deleting "<< bad_cells.size() << "bad cells" << endl;
  DeleteCells del;
  del.setGrid(m_Grid);
  del.setCellsToDelete(bad_cells);
  del();

  QMap<int,int> orgdir, curdir, voldir;
  EG_VTKDCC(vtkIntArray, cell_orgdir, m_Grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, m_Grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, m_Grid, "cell_voldir");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (cell_code->GetValue(id_cell) == bc) {
        orgdir[bc] = cell_orgdir->GetValue(id_cell);
        curdir[bc] = cell_curdir->GetValue(id_cell);
        voldir[bc] = cell_voldir->GetValue(id_cell);
        break;
      }
    }
    StitchHoles stitch(bc);
    stitch.setGrid(m_Grid);
    stitch();
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    int bc = cell_code->GetValue(id_cell);
    if (m_LayerAdjacentBoundaryCodes.contains(bc)) {
      cell_orgdir->SetValue(id_cell, orgdir[bc]);
      cell_curdir->SetValue(id_cell, curdir[bc]);
      cell_voldir->SetValue(id_cell, voldir[bc]);
    }
  }

  // check surface orientation
  {
    CorrectSurfaceOrientation surf_check;
    surf_check.setGrid(m_Grid);
    surf_check.setStart(0);
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (!m_LayerAdjacentBoundaryCodes.contains(cell_code->GetValue(id_cell))) {
        surf_check.setStart(id_cell);
        break;
      }
    }
    surf_check();
  }

  SwapTriangles swap;
  swap.setGrid(m_Grid);
  QSet<int> swap_codes = getAllBoundaryCodes(m_Grid);
  swap_codes -= m_LayerAdjacentBoundaryCodes;
  swap.setBoundaryCodes(swap_codes);
  swap.setVerboseOff();

  for (int iter = 0; iter < 5; ++iter) {
    cout << "correcting adjacent boundaries\n" << "  iteration: " << iter + 1 << endl;
    swap();
    smoothSurface();
    reduceSurface();
    swap();
  }
}
