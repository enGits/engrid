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
#include "pointfinder.h"

PointFinder::PointFinder()
{
  m_MinSize   = 1.0;
  m_MaxPoints = 100;
  m_SearchDistance = -1;
}

void PointFinder::setGrid(vtkUnstructuredGrid *grid)
{
  QVector<vec3_t> points(grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    grid->GetPoint(id_node, points[id_node].data());
  }
  setPoints(points);
}

void PointFinder::setPoints(const QVector<vec3_t> &points)
{
  m_Points = points;
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
    if (fabs(Dx1[0]) > fabs(Dx1[1]) && fabs(Dx1[0]) > fabs(Dx1[2])) {
      Dx1 = vec3_t(Dx1[0], Dx1[0], Dx1[0]);
    } else if (fabs(Dx1[1]) > fabs(Dx1[0]) && fabs(Dx1[1]) > fabs(Dx1[2])) {
      Dx1 = vec3_t(Dx1[1], Dx1[1], Dx1[1]);
    } else {
      Dx1 = vec3_t(Dx1[2], Dx1[2], Dx1[2]);
    }
    if (fabs(Dx2[0]) > fabs(Dx2[1]) && fabs(Dx2[0]) > fabs(Dx2[2])) {
      Dx2 = vec3_t(Dx2[0], Dx2[0], Dx2[0]);
    } else if (fabs(Dx2[1]) > fabs(Dx2[0]) && fabs(Dx2[1]) > fabs(Dx2[2])) {
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
}

int PointFinder::refine()
{
  int N = m_Octree.getNumCells();
  m_Octree.setSmoothTransitionOff();
  int N_new = 0;
  m_Octree.resetRefineMarks();
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    double dx = m_Octree.getDx(cell);
    if (!m_Octree.hasChildren(cell) && m_Buckets[cell].size() > m_MaxPoints  && dx > 2*m_MinSize) {
      m_Octree.markToRefine(cell);
      ++N_new;
    }
  }
  int N_mrk = 0;
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (m_Octree.markedForRefine(cell)) {
      ++N_mrk;
    }
  }
  N_new *= 8;
  m_Octree.refineAll();
  if (m_Octree.getNumCells() != N + N_new) {
    EG_BUG;
  }
  m_Buckets.insert(N, N_new, QList<int>());
  for (int cell = N; cell < m_Octree.getNumCells(); ++cell) {
    m_Timer << "  " << m_Octree.getNumCells() << " octree cells" << Timer::endl;
    if (m_Octree.getNumCells() != N + N_new) {
      EG_BUG;
    }
    int parent = m_Octree.getParent(cell);
    if (parent < 0) {
      EG_BUG;
    }
    foreach (int i_points, m_Buckets[parent]) {
      vec3_t xcell = m_Octree.getCellCentre(cell);
      vec3_t x = m_Points[i_points];
      bool append = true;
      if (m_SearchDistance > 0) {
        if      (x[0] < xcell[0] - m_SearchDistance || x[0] > xcell[0] + m_SearchDistance) append = false;
        else if (x[1] < xcell[1] - m_SearchDistance || x[1] > xcell[1] + m_SearchDistance) append = false;
        else if (x[2] < xcell[2] - m_SearchDistance || x[2] > xcell[2] + m_SearchDistance) append = false;
      } else {
        double Dx = m_Octree.getDx(cell);
        double Dy = m_Octree.getDy(cell);
        double Dz = m_Octree.getDz(cell);
        if      (x[0] < xcell[0] - Dx || x[0] > xcell[0] + Dx) append = false;
        else if (x[1] < xcell[1] - Dy || x[1] > xcell[1] + Dy) append = false;
        else if (x[2] < xcell[2] - Dz || x[2] > xcell[2] + Dz) append = false;
      }
      if (append) {
        m_Buckets[cell].append(i_points);
      }
    }
  }
  return N_new;
}

void PointFinder::getClosePoints(vec3_t x, QVector<int> &points)
{
  int cell = m_Octree.findCell(x);
  if (cell < 0) {
    EG_BUG;
  }
  if (m_Octree.hasChildren(cell)) {
    EG_BUG;
  }
  while (m_Buckets[cell].size() == 0 && m_Octree.getParent(cell) >= 0) {
    cell = m_Octree.getParent(cell);
  }
  points.resize(m_Buckets[cell].size());
  qCopy(m_Buckets[cell].begin(), m_Buckets[cell].end(), points.begin());

}

