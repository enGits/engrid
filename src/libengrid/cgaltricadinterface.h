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

  typedef QVector<Segment>::iterator                       SegmentIterator;
  typedef CGAL::AABB_segment_primitive<K, SegmentIterator> SegmentPrimitive;
  typedef CGAL::AABB_traits<K, SegmentPrimitive>           SegmentTraits;
  typedef CGAL::AABB_tree<SegmentTraits>                   SegmentTree;

private: // attributes

  QVector<Triangle> m_Triangles;
  TriangleTree      m_TriangleTree;
  QVector<Segment>  m_Segments;
  SegmentTree       m_SegmentTree;


public:

  CgalTriCadInterface(vtkUnstructuredGrid *grid);
  virtual vec3_t project(vec3_t, vec3_t, bool, bool) { notImplemented(); }
  virtual vec3_t snap(vec3_t x, bool correct_curvature = false);
  virtual vec3_t snapToEdge(vec3_t x);
  virtual vec3_t snapToCorner(vec3_t x);
  virtual double getRadius(vtkIdType id_node);


};

#endif // CGALTRICADINTERFACE_H
