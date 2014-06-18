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


#include "cadinterface.h"

class CgalTriCadInterface : public CadInterface
{

private: // types

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::FT FT;
  typedef K::Ray_3 Ray;
  typedef K::Line_3 Line;
  typedef K::Point_3 Point;
  typedef K::Triangle_3 Triangle;
  typedef QVector<Triangle>::iterator Iterator;
  typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


private: // attributes

  QVector<Triangle> m_Triangles;
  Tree              m_Tree;

public:

  CgalTriCadInterface(vtkUnstructuredGrid *grid);
  virtual vec3_t project(vec3_t, vec3_t, bool, bool) { notImplemented(); }
  virtual vec3_t  snap(vec3_t x, bool correct_curvature = false);


};

#endif // CGALTRICADINTERFACE_H
