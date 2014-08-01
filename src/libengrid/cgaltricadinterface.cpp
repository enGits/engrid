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

#include "cgaltricadinterface.h"
#include "meshpartition.h"

#include <CGAL/exceptions.h>

CgalTriCadInterface::CgalTriCadInterface(vtkUnstructuredGrid *grid)
{
  m_BGrid = vtkUnstructuredGrid::New();
  makeCopy(grid, m_BGrid);
  m_BPart.setGrid(m_BGrid, true);

  double feature_angle;

  EG_STOPDATE("2014-09-01");
  //getXmlSetting("engrid/surface/settings", "feature_angle", feature_angle);
  getSet("surface meshing", "feature angle", 20, feature_angle);
  feature_angle = GeometryTools::deg2rad(feature_angle);

  // build triangle tree
  {
    m_Triangles.clear();
    QVector<vtkIdType> tris;
    getAllCellsOfType(VTK_TRIANGLE, tris, m_BGrid);
    m_Triangles.fill(Triangle(), tris.size());
    m_Tri2Grid.fill(-1, tris.size());
    int i = 0;
    foreach (vtkIdType id_cell, tris) {
      EG_GET_CELL(id_cell, m_BGrid);
      vec3_t a, b, c;
      m_BGrid->GetPoint(pts[0], a.data());
      m_BGrid->GetPoint(pts[1], b.data());
      m_BGrid->GetPoint(pts[2], c.data());
      m_Triangles[i] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      m_Tri2Grid[i] = id_cell;
      ++i;
    }
    m_TriangleTree.rebuild(m_Triangles.begin(), m_Triangles.end());
    m_TriangleTree.accelerate_distance_queries();
  }

  // build edge tree
  {
    m_Segments.clear();
    QList<Segment> segs;
    EG_VTKDCC(vtkIntArray, cell_code, m_BGrid, "cell_code");
    EG_FORALL_NODES(id_node1, m_BGrid) {
      for (int i = 0; i < m_BPart.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_BPart.n2nGG(id_node1, i);
        if (id_node1 < id_node2) {
          QVector<vtkIdType> id_face;
          m_BPart.getEdgeFaces(id_node1, id_node2, id_face);
          if (id_face.size() == 2) {
            bool append_edge = false;
            if (cell_code->GetValue(id_face[0]) != cell_code->GetValue(id_face[1])) {
              append_edge = true;
            } else {
              vec3_t n1 = GeometryTools::cellNormal(m_BGrid, id_face[0]);
              vec3_t n2 = GeometryTools::cellNormal(m_BGrid, id_face[1]);
              double alpha = GeometryTools::angle(n1, n2);
              if (alpha >= feature_angle) {
                append_edge = true;
              }
            }
            if (append_edge) {
              vec3_t x1, x2;
              m_BGrid->GetPoint(id_node1, x1.data());
              m_BGrid->GetPoint(id_node2, x2.data());
              segs.append(Segment(Point(x1[0], x1[1], x1[2]), Point(x2[0], x2[1], x2[2])));
            }
          }
        }
      }
    }
    m_Segments.resize(segs.size());
    qCopy(segs.begin(), segs.end(), m_Segments.begin());
    m_SegmentTree.rebuild(m_Segments.begin(), m_Segments.end());
    m_SegmentTree.accelerate_distance_queries();
  }

  setName("CgalTriCadInterface");
  computeSurfaceCurvature();
  //m_ShootRayImplemented = true;
}

void CgalTriCadInterface::computeSurfaceCurvature()
{
  QVector<double> node_radius(m_BGrid->GetNumberOfPoints(), EG_LARGE_REAL);
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t x1;
    m_BGrid->GetPoint(id_node, x1.data());
    vec3_t n1 = m_BPart.globalNormal(id_node);
    for (int i = 0; i < m_BPart.n2nGSize(id_node); ++i) {
      vtkIdType id_neigh = m_BPart.n2nGG(id_node, i);
      vec3_t n2 = m_BPart.globalNormal(id_neigh);
      double scal_prod = max(-1.0, min(1.0, n1*n2));
      double alpha = max(1e-3, acos(scal_prod));
      if (alpha > 1e-3) {
        vec3_t x2;
        m_BGrid->GetPoint(id_neigh, x2.data());
        double a = (x1 - x2).abs();
        node_radius[id_node] = min(node_radius[id_node], a/alpha);
      }
    }
  }

  // compute weighted (distance) average of radii
  QVector<double> R_new(m_BGrid->GetNumberOfPoints(), 1e99);
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t x1;
    m_BGrid->GetPoint(id_node, x1.data());
    int N = m_BPart.n2nGSize(id_node);
    QVector<double> L(N);
    QVector<double> R(N, -1);
    double Lmax = 0;
    bool average = false;
    for (int i = 0; i < N; ++i) {
      vtkIdType id_neigh = m_BPart.n2nGG(id_node, i);
      vec3_t x2;
      m_BGrid->GetPoint(id_neigh, x2.data());
      L[i] = (x2 - x1).abs();
      if (node_radius[id_neigh] < 1e90 && L[i] > 0) {
        R[i] = node_radius[id_neigh];
        Lmax = max(Lmax, L[i]);
        average = true;
      }
    }
    if (average) {
      R_new[id_node] = 0;
      double total_weight = 0;
      for (int i = 0; i < N; ++i) {
        if (R[i] > 0) {
          R_new[id_node] += Lmax/L[i] * R[i];
          total_weight += Lmax/L[i];
          //R_new[id_node] += R[i];
          //total_weight += 1.0;
        }
      }
      R_new[id_node] /= total_weight;
    }
  }
  node_radius = R_new;

  m_Radius.fill(0, m_BGrid->GetNumberOfCells());
  EG_FORALL_CELLS(id_cell, m_BGrid) {
    EG_GET_CELL(id_cell, m_BGrid);
    for (int i = 0; i < num_pts; ++i) {
      m_Radius[id_cell] = max(m_Radius[id_cell], node_radius[pts[i]]);
    }
  }
}

CgalTriCadInterface::Ray CgalTriCadInterface::createRay(vec3_t x1, vec3_t v)
{
  Point p1(x1[0], x1[1], x1[2]);
  vec3_t x2 = x1 + v;
  Point p2(x2[0], x2[1], x2[2]);
  return Ray(p1, p2);
}

vec3_t CgalTriCadInterface::snap(vec3_t x, bool)
{
  vec3_t x_snap = x;
  try {
    Point p(x[0], x[1], x[2]);
    TrianglePointAndPrimitiveId result = m_TriangleTree.closest_point_and_primitive(p);
    Point cp = result.first;
    Triangle* T = result.second;
    int id = (T - m_Triangles.begin());
    vtkIdType id_face = m_Tri2Grid[id];
    m_LastNormal = GeometryTools::cellNormal(m_BGrid, id_face);
    m_LastNormal.normalise();
    m_LastRadius = m_Radius[id_face];
    x_snap = vec3_t(cp[0], cp[1], cp[2]);
  } catch (CGAL::Failure_exception) {
    x_snap = x;
  }
  return x_snap;
}

vec3_t CgalTriCadInterface::snapNode(vtkIdType id_node, vec3_t x, bool correct_curvature)
{
  vec3_t x_snap = x;
  try {
    x_snap = snap(x, correct_curvature);
    vec3_t n_node = m_FPart.globalNormal(id_node);
    if (n_node*m_LastNormal < 0) {
      double dist_min = EG_LARGE_REAL;
      Ray ray = createRay(x, n_node);
      std::list<Intersection> intersections;
      m_TriangleTree.all_intersections(ray, std::back_inserter(intersections));
      for (std::list<Intersection>::iterator i = intersections.begin(); i != intersections.end(); ++i) {
        int id = (i->second - m_Triangles.begin());
        vtkIdType id_face = m_Tri2Grid[id];
        vec3_t n = GeometryTools::cellNormal(m_BGrid, id_face);
        n.normalise();
        if (n*n_node > 0) {
          if (i->first.type() == typeid(Point)) {
            Point p = boost::get<Point>(i->first);
            vec3_t xs(p[0], p[1], p[2]);
            double dist = (x - xs).abs();
            if (dist < dist_min) {
              dist_min = dist;
              x_snap = xs;
              m_LastNormal = n;
              m_LastRadius = m_Radius[id_face];
            }
          }
        }
      }
    }
  } catch (CGAL::Failure_exception) {
    x_snap = x;
  }
  return x_snap;
}

vec3_t CgalTriCadInterface::snapToEdge(vec3_t x)
{
  vec3_t x_snap = x;
  try {
    Point p(x[0], x[1], x[2]);
    Point cp = m_SegmentTree.closest_point(p);
    x_snap = vec3_t(cp[0], cp[1], cp[2]);
  } catch (CGAL::Failure_exception) {
    x_snap = x;
  }
  return x_snap;
}

vec3_t CgalTriCadInterface::snapToCorner(vec3_t)
{
  notImplemented();
}

CgalTriCadInterface::HitType CgalTriCadInterface::shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r)
{
  notImplemented();
}

void CgalTriCadInterface::computeIntersections(vec3_t x, vec3_t v, QVector<QPair<vec3_t, vtkIdType> > &intersections)
{
  try {
    Ray ray = createRay(x, v);
    std::list<Intersection> inters;
    m_TriangleTree.all_intersections(ray, std::back_inserter(inters));
    int i_intersections = 0;
    for (std::list<Intersection>::iterator i = inters.begin(); i != inters.end(); ++i) {
      if (i->first.type() == typeid(Point)) {
        ++i_intersections;
      }
    }
    intersections.resize(i_intersections);
    i_intersections = 0;
    for (std::list<Intersection>::iterator i = inters.begin(); i != inters.end(); ++i) {
      int id = (i->second - m_Triangles.begin());
      if (i->first.type() == typeid(Point)) {
        intersections[i_intersections].second = m_Tri2Grid[id];
        Point p = boost::get<Point>(i->first);
        vec3_t xs(p[0], p[1], p[2]);
        intersections[i_intersections].first = xs;
        ++i_intersections;
      }
    }
  } catch (CGAL::Failure_exception) {
    intersections.clear();
  }
}
