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

CgalTriCadInterface::CgalTriCadInterface(vtkUnstructuredGrid *grid)
{
  double feature_angle;
  getSet("surface meshing", "feature angle", 20, feature_angle);
  feature_angle = GeometryTools::deg2rad(feature_angle);

  // build triangle tree
  {
    m_Triangles.clear();
    QVector<vtkIdType> tris;
    getAllCellsOfType(VTK_TRIANGLE, tris, grid);
    m_Triangles.fill(Triangle(), tris.size());
    int i = 0;
    foreach (vtkIdType id_cell, tris) {
      EG_GET_CELL(id_cell, grid);
      vec3_t a, b, c;
      grid->GetPoint(pts[0], a.data());
      grid->GetPoint(pts[1], b.data());
      grid->GetPoint(pts[2], c.data());
      m_Triangles[i] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      ++i;
    }
    m_TriangleTree.rebuild(m_Triangles.begin(), m_Triangles.end());
    m_TriangleTree.accelerate_distance_queries();
  }

  // build edge tree
  {
    m_Segments.clear();
    QList<Segment> segs;
    MeshPartition part(grid, true);
    EG_FORALL_NODES(id_node1, grid) {
      for (int i = 0; i < part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = part.n2nGG(id_node1, i);
        if (id_node1 < id_node2) {
          QVector<vtkIdType> id_face;
          part.getEdgeFaces(id_node1, id_node2, id_face);
          if (id_face.size() == 2) {
            vec3_t n1 = GeometryTools::cellNormal(grid, id_face[0]);
            vec3_t n2 = GeometryTools::cellNormal(grid, id_face[1]);
            double alpha = GeometryTools::angle(n1, n2);
            if (alpha >= feature_angle) {
              vec3_t x1, x2;
              grid->GetPoint(id_node1, x1.data());
              grid->GetPoint(id_node2, x2.data());
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
}

vec3_t CgalTriCadInterface::snap(vec3_t x, bool)
{
  Point p(x[0], x[1], x[2]);
  Point cp = m_TriangleTree.closest_point(p);
  vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
  return x_snap;
}

vec3_t CgalTriCadInterface::snapToEdge(vec3_t x)
{
  Point p(x[0], x[1], x[2]);
  Point cp = m_SegmentTree.closest_point(p);
  vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
  return x_snap;
}

vec3_t CgalTriCadInterface::snapToCorner(vec3_t)
{
  notImplemented();
}

double CgalTriCadInterface::getRadius(vtkIdType)
{
  /** @todo needs to be implemented
   *  http://doc.cgal.org/latest/Interpolation/index.html#secsurface
   */
  return EG_LARGE_REAL;
}
