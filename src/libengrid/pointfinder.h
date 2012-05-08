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
#ifndef POINTFINDER_H
#define POINTFINDER_H

#include "octree.h"
#include "triangle.h"
#include "timer.h"

#include <QVector>
#include <QList>

class PointFinder : public EgVtkObject
{

  Octree               m_Octree;
  QVector<vec3_t>      m_Points;
  double               m_MinSize;
  int                  m_MaxPoints;
  Timer                m_Timer;
  QVector<QList<int> > m_Buckets;
  int                  m_MinBucketSize;
  int                  m_MaxBucketSize;


private: // methods

  int refine();


public: // methods

  PointFinder();

  void setGrid(vtkUnstructuredGrid *grid);
  void setPoints(const QVector<vec3_t> &points);
  void setMaxNumPoints(int N) { m_MaxPoints = N; }
  void getClosePoints(vec3_t x, QVector<int> &points, double dist = 0);
  int  minBucketSize() { return m_MinBucketSize; }
  int  maxBucketSize() { return m_MaxBucketSize; }
  void writeOctreeMesh(QString file_name);

};


#endif // POINTFINDER_H
