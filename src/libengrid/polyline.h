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
#ifndef POLYLINE_H
#define POLYLINE_H

class PolyLine;

#include "curve.h"

#include <QVector>
#include <QFile>
#include <QTextStream>

class PolyLine : public Curve
{

protected: // attributes

  QVector<vec3_t> m_Points;
  double          m_Length;


protected: // methods

  template <typename Q_CONTAINER1, typename Q_CONTAINER2>
  int findClosest(vec3_t x, const Q_CONTAINER1& points, const Q_CONTAINER2& exclude);

  void getStencil(double l, int &i1, int &i2, double &w1, double &w2);

public:

  template <typename Q_CONTAINER>
  void fromPointsFile(QString points_file, const Q_CONTAINER& start_points);

  virtual vec3_t position(double l);
  virtual vec3_t normal(double l);
  virtual vec3_t intersection(vec3_t x, vec3_t n);

};

template <typename Q_CONTAINER>
void PolyLine::fromPointsFile(QString points_file, const Q_CONTAINER& start_points)
{
  // read raw point data from file
  QFile file(points_file);
  file.open(QIODevice::ReadOnly);
  QTextStream f(&file);
  QList<vec3_t> raw_points;
  while (!f.atEnd()) {
    QString line = f.readLine();
    QStringList parts = line.split(",");
    if (parts.size() >= 3) {
      vec3_t x;
      x[0] = parts[0].trimmed().toDouble();
      x[1] = parts[1].trimmed().toDouble();
      x[2] = parts[2].trimmed().toDouble();
      raw_points << x;
    }
  }

  // find start point of poly-line
  QList<int> exclude;
  int i_start = 0;
  double min_dist = EG_LARGE_REAL;
  foreach (vec3_t x, start_points) {
    int i_best = findClosest(x, raw_points, exclude);
    double dist = (x - raw_points[i_best]).abs();
    if (dist < min_dist) {
      min_dist = dist;
      i_start = i_best;
    }
  }

  // sort points in such a way that they form a curve ...
  m_Points.resize(raw_points.size());
  m_Points[0] = raw_points[i_start];
  exclude << i_start;
  for (int i = 1; i < raw_points.size(); ++i) {
    int i_next = findClosest(m_Points[i-1], raw_points, exclude);
    m_Points[i] = raw_points[i_next];
    exclude << i_next;
  }

  // compute the total length of the curve
  m_Length = 0;
  for (int i = 1; i < m_Points.size(); ++i) {
    m_Length += (m_Points[i] - m_Points[i-1]).abs();
  }
}

template <typename Q_CONTAINER1, typename Q_CONTAINER2>
int PolyLine::findClosest(vec3_t x, const Q_CONTAINER1& points, const Q_CONTAINER2& exclude)
{
  int i_best = -1;
  double min_dist = EG_LARGE_REAL;
  for (int i = 0; i < points.size(); ++i) {
    if (!exclude.contains(i)) {
      double dist = (points[i] - x).abs();
      if (dist < min_dist) {
        min_dist = dist;
        i_best = i;
      }
    }
  }
  return i_best;
}

#endif // POLYLINE_H
