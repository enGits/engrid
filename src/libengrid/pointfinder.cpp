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
#include "pointfinder.h"

PointFinder::PointFinder()
{
  m_MinSize   = 1.0;
  m_MaxPoints = 100;
}

void PointFinder::setGrid(vtkUnstructuredGrid *grid)
{
  QVector<vec3_t> points(grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    grid->GetPoint(id_node, points[id_node].data());
  }
  setPoints(points);
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
      {
        double Dx = m_Octree.getDx(cell);
        double Dy = m_Octree.getDy(cell);
        double Dz = m_Octree.getDz(cell);
        if      (x[0] < xcell[0] - 1.5*Dx || x[0] > xcell[0] + 1.5*Dx) append = false;
        else if (x[1] < xcell[1] - 1.5*Dy || x[1] > xcell[1] + 1.5*Dy) append = false;
        else if (x[2] < xcell[2] - 1.5*Dz || x[2] > xcell[2] + 1.5*Dz) append = false;
      }
      if (append) {
        m_Buckets[cell].append(i_points);
      }
    }
  }
  return N_new;
}

void PointFinder::getClosePoints(vec3_t x, QVector<int> &points, double dist)
{
  if (m_Octree.isInsideBounds(x)) {
    int cell = m_Octree.findCell(x);
    if (cell < 0) {
      EG_BUG;
    }
    if (m_Octree.hasChildren(cell)) {
      EG_BUG;
    }
    while ((m_Buckets[cell].size() == 0  || dist > m_Octree.getDx(cell)) && m_Octree.getParent(cell) >= 0) {
      cell = m_Octree.getParent(cell);
    }
    points.resize(m_Buckets[cell].size());
    qCopy(m_Buckets[cell].begin(), m_Buckets[cell].end(), points.begin());
  } else {
    points.resize(m_Points.size());
    for (int i = 0; i < points.size(); ++i) {
      points[i] = i;
    }
  }
}

void PointFinder::writeOctreeMesh(QString file_name)
{
  EG_VTKSP(vtkUnstructuredGrid, otg);
  m_Octree.toVtkGridHangingNodes(otg);
  writeGrid(otg, file_name);
}
