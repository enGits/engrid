// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "facefinder.h"
#include "triangle.h"

FaceFinder::FaceFinder()
{
  m_MinSize = 1.0;
  m_MaxFaces = 10;
  m_CollectedPointsIncrement = 10000;
  m_CollectedPoints.resize(10000);
  getSet("surface meshing", "use improved face search", false, m_UseImprovedFaceSearch);
}

void FaceFinder::setGrid(vtkUnstructuredGrid *grid)
{
  m_Grid = grid;
  m_Triangles.resize(m_Grid->GetNumberOfCells());
  m_Centres.resize(m_Grid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      EG_BUG;
    }
    m_Triangles[id_cell] = Triangle(m_Grid, id_cell);
    m_Centres[id_cell] = cellCentre(m_Grid, id_cell);
  }
  {
    double bounds[6];
    m_Grid->GetBounds(bounds);
    vec3_t x1(bounds[0], bounds[2], bounds[4]);
    vec3_t x2(bounds[1], bounds[3], bounds[5]);
    vec3_t xc = 0.5*(x1 + x2);
    vec3_t Dx = x2 - xc;
    double dx_max = max(Dx[0], max(Dx[1], Dx[2]));
    Dx = vec3_t(dx_max, dx_max, dx_max);
    x1 = xc - 2*Dx;
    x2 = xc + 2*Dx;
    m_Octree.setBounds(x1, x2);
    if (m_UseImprovedFaceSearch) {
      m_MinSize = 0.005*Dx[0];
    } else {
      m_MinSize = 0.0001*Dx[0];
    }
  }
  m_Faces.resize(1);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    m_Faces[0].append(id_cell);
  }  
  if (m_UseImprovedFaceSearch) {
    m_CritLengthOctreeCell.fill(EG_LARGE_REAL, 1);
    foreach (vtkIdType id_cell, m_Faces[0]) {
      m_CritLengthOctreeCell[0] = min(m_CritLengthOctreeCell[0], calcCritLength(id_cell));
    }
    calcCritLengthForAllNodes();
  } else {
    m_CritLengthOctreeCell.fill(0, 1);
    foreach (vtkIdType id_cell, m_Faces[0]) {
      m_CritLengthOctreeCell[0] = max(m_CritLengthOctreeCell[0], calcCritLength(id_cell));
    }
  }

  int N;
  do {
    N = refine();
  } while (N > 0);

  int max_num_faces = 0;
  double ave_num_faces;
  int num_buckets = 0;
  for (int i_cell = 0; i_cell < m_Octree.getNumCells(); ++i_cell) {
    if (!m_Octree.hasChildren(i_cell)) {
      int cell = i_cell;
      while (m_Faces[cell].size() == 0) {
        cell = m_Octree.getParent(cell);
        if (cell == -1) {
          EG_BUG;
        }
      }
      max_num_faces = max(m_Faces[cell].size(), max_num_faces);
      ave_num_faces += double(m_Faces[cell].size());
      ++num_buckets;
    }
  }

  cout << "refining octree for FaceFinder ..." << endl;
  cout << "  minimal cell size: " << m_MinSize << endl;
  cout << "  maximal number of faces: " << max_num_faces << endl;
  cout << "  average number of faces: " << ave_num_faces/num_buckets << endl;
  cout << "  number of buckets: " << num_buckets << endl;
  //cout << "  largest critical length: " << m_CritLength[0] << endl;

  //EG_VTKSP(vtkUnstructuredGrid, otg);
  //m_Octree.toVtkGrid(otg);
  //writeGrid(otg, "facefinder");
}

double FaceFinder::calcCritLength(vtkIdType id_cell)
{
  if (m_UseImprovedFaceSearch) {
    return max(m_MinSize, 10*m_Triangles[id_cell].smallestLength());
  }
  QVector<vec3_t> x;
  vtkIdType N_pts, *pts;
  m_Grid->GetCellPoints(id_cell, N_pts, pts);
  x.resize(N_pts + 1);
  for (int i = 0; i < N_pts; ++i) {
    m_Grid->GetPoint(pts[i], x[i].data());
  }
  x[N_pts] = x[0];
  double L = 0;
  for (int i = 0; i < N_pts; ++i) {
    L = max(L, (x[i]-x[i+1]).abs());
  }
  return L;
}

void FaceFinder::calcCritLengthForAllNodes()
{
  m_CritLengthNode.fill(EG_LARGE_REAL, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_face = 0; id_face < m_Grid->GetNumberOfCells(); ++id_face) {
    double l_crit = calcCritLength(id_face);
    vtkIdType num_pts, *pts;
    m_Grid->GetCellPoints(id_face, num_pts, pts);
    for (int i = 0; i < num_pts; ++i) {
      m_CritLengthNode[pts[i]] = min(m_CritLengthNode[pts[i]], l_crit);
    }
  }
}

void FaceFinder::getPointsOfFace(vtkIdType id_face)
{
  m_NumCollectedPoints = 3;
  Triangle T = m_Triangles[id_face];
  m_CollectedPoints[0] = T.a();
  m_CollectedPoints[1] = T.b();
  m_CollectedPoints[2] = T.c();
  QVector<vec3_t> x_tri(4);
  QVector<vtkIdType> id_tri(4);
  x_tri[0] = T.a();
  x_tri[1] = T.b();
  x_tri[2] = T.c();
  x_tri[3] = T.a();
  id_tri[0] = T.idA();
  id_tri[1] = T.idB();
  id_tri[2] = T.idC();
  id_tri[3] = T.idA();
  for (int i = 0; i < 3; ++i) {
    vec3_t x1 = x_tri[i];
    vec3_t x2 = x_tri[i+1];
    double l  = (x2 - x1).abs();
    double h = 0.5*(m_CritLengthNode[id_tri[i]] + m_CritLengthNode[id_tri[i+1]]);
    int N = int(l/h) + 1;
    vec3_t dx = (1.0/N)*(x2 - x1);
    vec3_t x = x1 + dx;
    for (int j = 1; j < N; ++j) {
      if (m_NumCollectedPoints >= m_CollectedPoints.size()) {
        m_CollectedPoints.resize(m_CollectedPoints.size() + m_CollectedPointsIncrement);
      }
      m_CollectedPoints[m_NumCollectedPoints] = x;
      ++m_NumCollectedPoints;
      x += dx;
    }
  }
}

int FaceFinder::refine()
{
  int N = m_Octree.getNumCells();
  m_Octree.setSmoothTransitionOff();
  int N_new = 0;
  m_Octree.resetRefineMarks();
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    double dx = m_Octree.getDx(cell);
    if (!m_Octree.hasChildren(cell)  && m_Faces[cell].size() > m_MaxFaces  && dx > 2*m_MinSize && dx > 2*m_CritLengthOctreeCell[cell]) {
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
  m_Faces.insert(N, N_new, QList<vtkIdType>());
  if (m_UseImprovedFaceSearch) {
    m_CritLengthOctreeCell.insert(N, N_new, 1e99);
  } else {
    m_CritLengthOctreeCell.insert(N, N_new, 0);
  }
  int N_min = m_Grid->GetNumberOfCells();
  int N_max = 0;
  int N_ave = 0;
  for (int cell = N; cell < m_Octree.getNumCells(); ++cell) {    
    m_Timer << "  " << m_Octree.getNumCells() << " octree cells" << Timer::endl;
    if (m_Octree.getNumCells() != N + N_new) {
      EG_BUG;
    }
    int parent = m_Octree.getParent(cell);
    if (parent < 0) {
      EG_BUG;
    }
    foreach (vtkIdType id_cell, m_Faces[parent]) {
      if (m_UseImprovedFaceSearch) {
        getPointsOfFace(id_cell);
        for (int i = 0; i < m_NumCollectedPoints; ++i) {
          if (m_Octree.isInsideCell(cell, m_CollectedPoints[i], 2.0)) {
            m_Faces[cell].append(id_cell);
            m_CritLengthOctreeCell[cell] = min(m_CritLengthOctreeCell[cell], calcCritLength(id_cell));
            break;
          }
        }
      } else {
        for (int node = 0; node < 8; ++node) {
          vec3_t x = m_Octree.getNodePosition(cell, node);
          double d;
          d = (x - m_Centres[id_cell]).abs();
          if (d < m_Octree.getDx(cell)) {
            m_Faces[cell].append(id_cell);
            m_CritLengthOctreeCell[cell] = max(m_CritLengthOctreeCell[cell], calcCritLength(id_cell));
            break;
          }
        }
      }
    }
    N_min = min(N_min, m_Faces[cell].size());
    N_max = max(N_max, m_Faces[cell].size());
    N_ave += m_Faces[cell].size();
  }
  int N_cells = 0;
  int N_faces = 0;
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (!m_Octree.hasChildren(cell)) {
      ++N_cells;
      N_faces += m_Faces[cell].size();
    }
  }
  if (N_faces == 0) {
    EG_BUG;
  }
  return N_new;
}

void FaceFinder::getCloseFaces(vec3_t x, QVector<vtkIdType> &faces)
{
  if (m_Octree.isInsideBounds(x)) {
    int cell = m_Octree.findCell(x);
    if (cell < 0) {
      EG_BUG;
    }
    if (m_Octree.hasChildren(cell)) {
      EG_BUG;
    }
    while (m_Faces[cell].size() == 0 && m_Octree.getParent(cell) >= 0) {
      cell = m_Octree.getParent(cell);
    }
    faces.resize(m_Faces[cell].size());
    qCopy(m_Faces[cell].begin(), m_Faces[cell].end(), faces.begin());
  } else {
    faces.clear();
  }
}

vtkIdType FaceFinder::getClosestFace(vec3_t x, double &L_min)
{
  QVector<vtkIdType> faces;
  getCloseFaces(x, faces);
  vtkIdType id_close = -1;
  foreach (vtkIdType id_face, faces) {
    vec3_t xi, ri;
    int side;
    double L;
    m_Triangles[id_face].snapOntoTriangle(x, xi, ri, L, side, true);
    if (L < L_min) {
      L_min = L;
      id_close = id_face;
    }
  }
  return id_close;
}

