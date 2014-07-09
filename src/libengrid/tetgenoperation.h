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
// 
#ifndef TETGENOPERATION_H
#define TETGENOPERATION_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

#include "operation.h"
#include "tetgen.h"
#include "edgelengthsourcemanager.h"

class TetGenOperation;

class TetGenOperation : public Operation
{

protected: // data types

  struct segment_t
  {
    int node1, node2;
  };

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


protected: // attributes

  double m_MinimalEdgeLength;
  double m_MaximalEdgeLength;
  double m_GrowthFactor;
  double m_NodesPerQuarterCircle;
  double m_2dFeatureResolution;
  double m_3dFeatureResolution;
  double m_FeatureAngle;
  int    m_OrgDir;
  int    m_CurDir;
  int    m_VolDir;

  EdgeLengthSourceManager m_ELSManager;

  QString m_TetGenPath;


protected: // methods

  void copyToTetGen(tetgenio &tgio, vtkUnstructuredGrid *alt_grid = NULL);
  void copyFromTetGen(tetgenio &tgio);

  void tetgen(QString flags, vtkUnstructuredGrid *background_grid = NULL);
  void readSettings();

  QString qualityText();


public:

  TetGenOperation();

};

#endif // TETGENOPERATION_H
