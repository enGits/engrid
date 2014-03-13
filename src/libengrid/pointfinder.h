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

  template <typename C>
  void setPoints(const C &points);

  void setGrid(vtkUnstructuredGrid *grid);
  void setMaxNumPoints(int N) { m_MaxPoints = N; }
  void getClosePoints(vec3_t x, QVector<int> &points, double dist = 0);
  int  minBucketSize() { return m_MinBucketSize; }
  int  maxBucketSize() { return m_MaxBucketSize; }
  void writeOctreeMesh(QString file_name);

};



template <typename C>
void PointFinder::setPoints(const C &points)
{
  m_Points.resize(points.size());
  qCopy(points.begin(), points.end(), m_Points.begin());
  {
    vec3_t x1(1e99, 1e99, 1e99);
    vec3_t x2(-1e99, -1e99, -1e99);
    foreach (vec3_t x, m_Points) {
      x1[0] = min(x[0], x1[0]);
      x1[1] = min(x[1], x1[1]);
      x1[2] = min(x[2], x1[2]);
      x2[0] = max(x[0], x2[0]);
      x2[1] = max(x[1], x2[1]);
      x2[2] = max(x[2], x2[2]);
    }
    vec3_t xc = 0.5*(x1 + x2);
    vec3_t Dx1 = xc - x1;
    vec3_t Dx2 = x2 - xc;
    if (fabs(Dx1[0]) >= fabs(Dx1[1]) && fabs(Dx1[0]) >= fabs(Dx1[2])) {
      Dx1 = vec3_t(Dx1[0], Dx1[0], Dx1[0]);
    } else if (fabs(Dx1[1]) >= fabs(Dx1[0]) && fabs(Dx1[1]) >= fabs(Dx1[2])) {
      Dx1 = vec3_t(Dx1[1], Dx1[1], Dx1[1]);
    } else {
      Dx1 = vec3_t(Dx1[2], Dx1[2], Dx1[2]);
    }
    if (fabs(Dx2[0]) >= fabs(Dx2[1]) && fabs(Dx2[0]) >= fabs(Dx2[2])) {
      Dx2 = vec3_t(Dx2[0], Dx2[0], Dx2[0]);
    } else if (fabs(Dx2[1]) >= fabs(Dx2[0]) && fabs(Dx2[1]) >= fabs(Dx2[2])) {
      Dx2 = vec3_t(Dx2[1], Dx2[1], Dx2[1]);
    } else {
      Dx2 = vec3_t(Dx2[2], Dx2[2], Dx2[2]);
    }
    vec3_t Dx = Dx1;
    if (Dx2.abs() > Dx1.abs()) {
      Dx = Dx2;
    }
    x1 = xc - 2*Dx;
    x2 = xc + 2*Dx;
    m_Octree.setBounds(x1, x2);
    m_MinSize = 0.0001*Dx[0];
  }
  m_Buckets.resize(1);
  for (int i_points = 0; i_points < points.size(); ++i_points) {
    m_Buckets[0].append(i_points);
  }
  int N;
  do {
    N = refine();
  } while (N > 0);
  m_MaxBucketSize = 0;
  m_MinBucketSize = points.size();
  for (int i = 0; i < m_Buckets.size(); ++i) {
    m_MinBucketSize = min(m_MinBucketSize, m_Buckets[i].size());
    m_MaxBucketSize = max(m_MaxBucketSize, m_Buckets[i].size());
  }
}



#endif // POINTFINDER_H
