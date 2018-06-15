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

#ifndef CGALTRICADINTERFACE_H
#define CGALTRICADINTERFACE_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>


#include "cadinterface.h"

class CgalTriCadInterface : public CadInterface
{

private: // types

  typedef CGAL::Simple_cartesian<double> K;

  typedef K::FT         FT;
  typedef K::Ray_3      Ray;
  typedef K::Line_3     Line;
  typedef K::Point_3    Point;
  typedef K::Segment_3  Segment;
  typedef K::Triangle_3 Triangle;

  typedef QVector<Triangle>::iterator                       TriangleIterator;
  typedef CGAL::AABB_triangle_primitive<K,TriangleIterator> TrianglePrimitive;
  typedef CGAL::AABB_traits<K, TrianglePrimitive>           TriangleTraits;
  typedef CGAL::AABB_tree<TriangleTraits>                   TriangleTree;
  typedef TriangleTree::Point_and_primitive_id              TrianglePointAndPrimitiveId;

  typedef QVector<Segment>::iterator                       SegmentIterator;
  typedef CGAL::AABB_segment_primitive<K, SegmentIterator> SegmentPrimitive;
  typedef CGAL::AABB_traits<K, SegmentPrimitive>           SegmentTraits;
  typedef CGAL::AABB_tree<SegmentTraits>                   SegmentTree;

  //typedef boost::optional<TriangleTree::Intersection_and_primitive_id<Ray>::Type> Intersection;
  typedef TriangleTree::Intersection_and_primitive_id<Ray>::Type Intersection;

private: // attributes

  QVector<Triangle>    m_Triangles;
  TriangleTree         m_TriangleTree;
  QVector<Segment>     m_Segments;
  SegmentTree          m_SegmentTree;
  vtkUnstructuredGrid *m_BGrid;
  MeshPartition        m_BPart;
  QVector<vtkIdType>   m_Tri2Grid;
  QVector<double>      m_Radius;      ///< Surface radius for mesh resolution.
  vtkIdType            m_LastFaceId;


private: // methods

  void computeSurfaceCurvature();
  Ray createRay(vec3_t x1, vec3_t v);


public:

  CgalTriCadInterface(vtkUnstructuredGrid *grid);
  virtual HitType shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r);
  virtual vec3_t snap(vec3_t x, bool correct_curvature = false);
  virtual vec3_t snapWithNormal(vec3_t x, vec3_t n, bool correct_curvature = false);
  virtual vec3_t snapNode(vtkIdType id_node, vec3_t x, bool correct_curvature);
  virtual vec3_t snapToEdge(vec3_t x);
  virtual vec3_t snapToCorner(vec3_t x);
  virtual void computeIntersections(vec3_t x, vec3_t v, QVector<QPair<vec3_t, vtkIdType> > &intersections);

  vtkIdType getLastFaceId() { return m_LastFaceId; }

};

#endif // CGALTRICADINTERFACE_H
