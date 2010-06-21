//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "surfaceprojection.h"

#include <math.h>

long int SurfaceProjection::Nfull = 0;
long int SurfaceProjection::Nhalf = 0;
vtkIdType SurfaceProjection::m_LastPindex = 0;

SurfaceProjection::SurfaceProjection(int bc) : SurfaceAlgorithm()
{
  m_BGrid = vtkUnstructuredGrid::New();
  this->setGrid(m_BGrid);
  m_BC = bc;
  getSet("surface meshing", "correct curvature (experimental)", false, m_correctCurvature);
  getSet("surface meshing", "number ofneighbour searches", 2, m_NumNeighSearches);
  m_CritDistance = 0.1;
}

SurfaceProjection::~SurfaceProjection()
{
  m_BGrid->Delete();
}

void SurfaceProjection::setForegroundGrid(vtkUnstructuredGrid *grid)
{
  m_FGrid = grid;
}

void SurfaceProjection::searchNewTriangle(vec3_t xp, vtkIdType &id_tri, vec3_t &x_proj, vec3_t &r_proj, bool &on_triangle)
{
  x_proj = vec3_t(1e99, 1e99, 1e99);
  r_proj = vec3_t(0, 0, 0);
  vec3_t xi, ri;
  on_triangle = false;
  double d_min = 1e99;
  bool tri_found = false;
  if (id_tri >= 0) {
    QSet<vtkIdType> tris;
    tris.insert(id_tri);
    for (int i_search = 0; i_search < m_NumNeighSearches; ++i_search) {
    //for (int i_search = 0; i_search < 100; ++i_search) {
      QSet<vtkIdType> candidates;
      foreach (int i_tri, tris) {
        for (int k = 0; k < m_BPart.c2cGSize(i_tri); ++k) {
          int new_cand = m_BPart.c2cGG(i_tri,k);
          if (!tris.contains(new_cand) && new_cand >= 0) {
            candidates.insert(m_BPart.c2cGG(i_tri,k));
          }
        }
      }
      foreach (int i_cand, candidates) {
        if (i_cand < 0) EG_BUG;
        if (i_cand > m_BGrid->GetNumberOfCells()) EG_BUG;
        int side;
        double d;
        if(m_Triangles[i_cand].projectOnTriangle(xp, xi, ri, d, side, m_RestrictToTriangle)) {
          x_proj = xi;
          r_proj = ri;
          on_triangle = true;
          id_tri = i_cand;
          ++Nhalf;
          return;
        } else if (d < m_CritDistance*m_Triangles[i_cand].m_smallest_length) {
          if (d < d_min) {
            x_proj = xi;
            r_proj = ri;
            id_tri = i_cand;
            tri_found = true;
          }
        }
      }
      tris += candidates;
      if (tri_found) {
        return;
      }
    }
  }
  if (tri_found) {
    ++Nhalf;
  } else {
    // full search :-(
    ++Nfull;
    bool x_proj_set = false;
    d_min = 1e99;
    for (int i_triangles = 0; i_triangles < m_Triangles.size(); ++i_triangles) {
      Triangle T = m_Triangles[i_triangles];
      double d;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, m_RestrictToTriangle);
      if (d >= 1e99) {
        EG_BUG;
      }
      if (d < d_min) {
        x_proj = xi;
        x_proj_set = true;
        r_proj = ri;
        d_min = d;
        id_tri = i_triangles;
        on_triangle = intersects;
      }
    }
    if (!x_proj_set) { // should never happen
      checkVector(xp);
      qWarning() << "No projection found for point xp=" << xp[0] << xp[1] << xp[2] << endl;
      writeGrid(GuiMainWindow::pointer()->getGrid(), "griddump");
      EG_BUG;
    }
  }
}

void SurfaceProjection::updateBackgroundGridInfo()
{
  getAllCells(m_Cells, m_BGrid);
  getNodesFromCells(m_Cells, m_Nodes, m_BGrid);
  setBoundaryCodes(GuiMainWindow::pointer()->getAllBoundaryCodes());
  setAllCells();
  readVMD();
  UpdatePotentialSnapPoints(true, false);
  l2l_t  c2c   = getPartC2C();
  g2l_t _cells = getPartLocalCells();
  QVector<int> m_LNodes(m_Nodes.size());
  for (int i = 0; i < m_LNodes.size(); ++i) {
    m_LNodes[i] = i;
  }
  createNodeToNode(m_Cells, m_Nodes, m_LNodes, m_N2N, m_BGrid);

  // create m_Triangles
  m_Triangles.resize(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    m_Triangles[id_cell] = Triangle(m_BGrid, id_cell);
    for (int i = 0; i < 3; i++) {
      int i_cell = _cells[id_cell];
      if (c2c[i_cell][i] < 0) {
        m_Triangles[id_cell].m_has_neighbour[i] = false;
      } else {
        m_Triangles[id_cell].m_has_neighbour[i] = true;
      }
    }
  }

  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0, 0, 0);
  }

  foreach(Triangle T, m_Triangles) {
    double angle_a = GeometryTools::angle(m_BGrid, T.m_id_c, T.m_id_a, T.m_id_b);
    double angle_b = GeometryTools::angle(m_BGrid, T.m_id_a, T.m_id_b, T.m_id_c);
    double angle_c = GeometryTools::angle(m_BGrid, T.m_id_b, T.m_id_c, T.m_id_a);
    if (isnan(angle_a) || isinf(angle_a)) EG_BUG;
    if (isnan(angle_b) || isinf(angle_b)) EG_BUG;
    if (isnan(angle_c) || isinf(angle_c)) EG_BUG;
    if (!checkVector(T.m_g3)) {
      qWarning() << "T.m_g3=" << T.m_g3;
      EG_BUG;
    }
    m_NodeNormals[T.m_id_a] += angle_a * T.m_g3;
    m_NodeNormals[T.m_id_b] += angle_b * T.m_g3;
    m_NodeNormals[T.m_id_c] += angle_c * T.m_g3;
    if (!checkVector(m_NodeNormals[T.m_id_a])) EG_BUG;
    if (!checkVector(m_NodeNormals[T.m_id_b])) EG_BUG;
    if (!checkVector(m_NodeNormals[T.m_id_c])) EG_BUG;
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node].normalise();
  }
  m_BPart.setGrid(m_BGrid);
  m_BPart.setAllCells();
}


vtkIdType SurfaceProjection::getProjTriangle(vtkIdType id_node)
{
  EG_VTKDCN(vtkLongArray_t, pi, m_FGrid, "node_pindex");
  vtkIdType pindex = pi->GetValue(id_node);
  vtkIdType proj_triangle = -1;
  if (pindex < 0) {
    pindex = m_LastPindex;
    pi->SetValue(id_node, pindex);
    ++m_LastPindex;
    m_Pindex[pindex] = -1;
  } else {
    proj_triangle = m_Pindex[pindex];
  }
  return proj_triangle;
}

void SurfaceProjection::setProjTriangle(vtkIdType id_node, vtkIdType proj_triangle)
{
  EG_VTKDCN(vtkLongArray_t, pi, m_FGrid, "node_pindex");
  vtkIdType pindex = pi->GetValue(id_node);
  if (pindex < 0) {
    pindex = m_LastPindex;
    pi->SetValue(id_node, pindex);
    ++m_LastPindex;
  }
  m_Pindex[pindex] = proj_triangle;
}


vec3_t SurfaceProjection::project(vec3_t xp, vtkIdType id_node)
{
  checkVector(xp);
  if (id_node < 0) {
    EG_BUG;
  }

  vec3_t x_proj(1e99, 1e99, 1e99);
  vec3_t r_proj(0, 0, 0);
  bool x_proj_set = false;

  // initilizing booleans
  bool on_triangle = false;
  bool need_search = false;

  vtkIdType proj_triangle = getProjTriangle(id_node);

  if (proj_triangle == -1) { //if there is no known triangle on which to project
    need_search = true;
  } else { //if there is a known triangle on which to project
    if (proj_triangle >= m_Triangles.size()) {
      EG_BUG;
    }
    Triangle T = m_Triangles[proj_triangle];
    vec3_t xi, ri;
    double d;
    int side;
    bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, m_RestrictToTriangle);
    if (!intersects || (d > m_CritDistance*T.m_smallest_length)) {
      need_search = true;
    } else {
      x_proj = xi;
      x_proj_set = true;
      if (x_proj[0] > 1e98) { // should never happen
        EG_BUG;
      }
      r_proj = ri;
      on_triangle = intersects;
    }
  }
  //need_search = true;
  if (need_search) {
    searchNewTriangle(xp, proj_triangle, x_proj, r_proj, on_triangle);
    setProjTriangle(id_node, proj_triangle);
  }
  if (m_correctCurvature) {
    x_proj = correctCurvature(proj_triangle, xp);
  }
  return x_proj;
}

vec3_t SurfaceProjection::projectFree(vec3_t x, vtkIdType id_node)
{
  m_RestrictToTriangle = false;
  return project(x, id_node);
}

vec3_t SurfaceProjection::projectRestricted(vec3_t x, vtkIdType id_node)
{
  m_RestrictToTriangle = true;
  return project(x, id_node);
}

vec3_t SurfaceProjection::correctCurvature(int, vec3_t g_M)
{
  return g_M;
}
