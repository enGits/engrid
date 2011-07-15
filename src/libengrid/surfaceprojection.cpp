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
  double d_min      = 1e99;
  bool   x_proj_set = false;
  on_triangle = false;
  QVector<vtkIdType> close_faces;
  m_FaceFinder.getCloseFaces(xp, close_faces);
  foreach (vtkIdType id_triangle, close_faces) {
    Triangle T = m_Triangles[id_triangle];
    double d;
    int side;
    bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
    if (d >= 1e99) {
      EG_BUG;
    }
    if (d < d_min) {
      x_proj = xi;
      x_proj_set = true;
      r_proj = ri;
      d_min = d;
      id_tri = id_triangle;
      on_triangle = intersects;
    }
  }
  if (!x_proj_set) { // should never happen
    for (vtkIdType id_triangle = 0; id_triangle < m_BGrid->GetNumberOfCells(); ++id_triangle) {
      Triangle T = m_Triangles[id_triangle];
      double d;
      int side;
      bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, true);
      if (d >= 1e99) {
        EG_BUG;
      }
      if (d < 0) {
        EG_BUG;
      }
      if (d < d_min) {
        x_proj = xi;
        x_proj_set = true;
        r_proj = ri;
        d_min = d;
        id_tri = id_triangle;
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
        m_Triangles[id_cell].setNeighbourFalse(i);
      } else {
        m_Triangles[id_cell].setNeighbourTrue(i);
      }
    }
  }

  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0, 0, 0);
  }

  foreach(Triangle T, m_Triangles) {
    double angle_a = GeometryTools::angle(m_BGrid, T.idC(), T.idA(), T.idB());
    double angle_b = GeometryTools::angle(m_BGrid, T.idA(), T.idB(), T.idC());
    double angle_c = GeometryTools::angle(m_BGrid, T.idB(), T.idC(), T.idA());
    if (isnan(angle_a) || isinf(angle_a)) EG_BUG;
    if (isnan(angle_b) || isinf(angle_b)) EG_BUG;
    if (isnan(angle_c) || isinf(angle_c)) EG_BUG;
    if (!checkVector(T.g3())) {
      qWarning() << "T.g3 = " << T.g3();
      EG_BUG;
    }
    m_NodeNormals[T.idA()] += angle_a * T.g3();
    m_NodeNormals[T.idB()] += angle_b * T.g3();
    m_NodeNormals[T.idC()] += angle_c * T.g3();
    if (!checkVector(m_NodeNormals[T.idA()])) EG_BUG;
    if (!checkVector(m_NodeNormals[T.idB()])) EG_BUG;
    if (!checkVector(m_NodeNormals[T.idC()])) EG_BUG;
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node].normalise();
  }
  for (int i = 0; i < m_Triangles.size(); ++i) {
    Triangle &T = m_Triangles[i];
    T.setNormals(m_NodeNormals[T.idA()], m_NodeNormals[T.idB()], m_NodeNormals[T.idC()]);
  }

  m_BPart.setGrid(m_BGrid);
  m_BPart.setAllCells();
  computeSurfaceCurvature();
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


vec3_t SurfaceProjection::project(vec3_t xp, vtkIdType id_node,  bool correct_curvature)
{
  if (!checkVector(xp)) {
    qWarning() << "No projection found for point, id_node=" << id_node << ", xp=" << xp[0] << xp[1] << xp[2] << endl;
    writeGrid(GuiMainWindow::pointer()->getGrid(), "griddump");
    EG_BUG;
  }

  vec3_t x_proj(1e99, 1e99, 1e99);
  vec3_t r_proj(0, 0, 0);
  bool x_proj_set = false;

  // initilizing booleans
  bool on_triangle = false;

  vtkIdType proj_triangle = -1;
  if (id_node != -1) {
    proj_triangle = getProjTriangle(id_node);
  }

  if (proj_triangle == -1) {
    searchNewTriangle(xp, proj_triangle, x_proj, r_proj, on_triangle);
    if (id_node != -1) {
      setProjTriangle(id_node, proj_triangle);
    }
  }
  if (proj_triangle >= m_Triangles.size()) {
    EG_BUG;
  }
  Triangle T = m_Triangles[proj_triangle];
  vec3_t xi, ri;
  double d;
  int side;
  bool intersects = T.projectOnTriangle(xp, xi, ri, d, side, m_RestrictToTriangle);
  if (!intersects || (d > m_CritDistance*T.smallestLength())) {
    searchNewTriangle(xp, proj_triangle, x_proj, r_proj, on_triangle);
    T = m_Triangles[proj_triangle];
    T.projectOnTriangle(xp, xi, ri, d, side, m_RestrictToTriangle);
    if (id_node != -1) {
      setProjTriangle(id_node, proj_triangle);
    }
  }
  x_proj = xi;
  x_proj_set = true;
  if (x_proj[0] > 1e98) { // should never happen
    EG_BUG;
  }
  r_proj = ri;
  on_triangle = intersects;
  if (!checkVector(x_proj)) {
    EG_BUG;
  }
  if (correct_curvature) {
    vec3_t x_corr = correctCurvature(proj_triangle, x_proj);
    if (checkVector(x_corr)) {
      x_proj = x_corr;
    }
  }
  m_LastProjTriangle = proj_triangle;
  return x_proj;
}

vec3_t SurfaceProjection::projectFree(vec3_t x, vtkIdType id_node, bool correct_curvature)
{
  m_RestrictToTriangle = false;
  return project(x, id_node, correct_curvature);
}

vec3_t SurfaceProjection::projectRestricted(vec3_t x, vtkIdType id_node, bool correct_curvature)
{
  m_RestrictToTriangle = true;
  return project(x, id_node, correct_curvature);
}

vec3_t SurfaceProjection::correctCurvature(vtkIdType proj_triangle, vec3_t x)
{
  vec3_t x_corr = x;
  if (proj_triangle != -1) {
    Triangle &T = m_Triangles[proj_triangle];
    vec3_t rx = T.global3DToLocal3D(x);
    double w1 = 1.0 - rx[0] - rx[1];
    double w2 = rx[0];
    double w3 = rx[1];
    double k1 = GeometryTools::intersection(x, T.g3(), T.a(), T.nA());
    double k2 = GeometryTools::intersection(x, T.g3(), T.b(), T.nB());
    double k3 = GeometryTools::intersection(x, T.g3(), T.c(), T.nC());
    double S = 0.5;
    double k = w1*k1 + w2*k2 + w3*k3;
    k -= rx[2];
    x_corr = x + S*k*T.g3();
    if (!checkVector(x_corr)) {
      x_corr = x;
    }
  }
  return x_corr;
}

void SurfaceProjection::computeSurfaceCurvature()
{
  /*
  m_Radius.fill(1e99, m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_BGrid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vec3_t> x(N_pts+1);
    QVector<vec3_t> n(N_pts+1);
    for (int i = 0; i < N_pts; ++i) {
      m_BGrid->GetPoint(pts[i], x[i].data());
      n[i] = m_NodeNormals[pts[i]];
    }
    x[N_pts] = x[0];
    n[N_pts] = n[0];
    for (int i = 0; i < N_pts; ++i) {
      double alpha = max(1e-3,acos(n[i]*n[i+1]));
      if (alpha > 1e-3) {
        double a = (x[i]-x[i+1]).abs();
        m_Radius[id_cell] = min(m_Radius[id_cell], a/alpha);
      }
    }
  }
  */
  m_Radius.fill(1e99, m_BGrid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t x1;
    m_BGrid->GetPoint(id_node, x1.data());
    for (int i = 0; i < m_BPart.n2nGSize(id_node); ++i) {
      vtkIdType id_neigh = m_BPart.n2nGG(id_node, i);
      double scal_prod = max(-1.0, min(1.0, m_NodeNormals[id_node]*m_NodeNormals[id_neigh]));
      double alpha = max(1e-3, acos(scal_prod));
      if (alpha > 1e-3) {
        vec3_t x2;
        m_BGrid->GetPoint(id_neigh, x2.data());
        double a = (x1 - x2).abs();
        m_Radius[id_node] = min(m_Radius[id_node], a/alpha);
        if (m_Radius[id_node] < 1e-3) {
          cout << "Science goes Boink!" << endl;
        }
      }
    }
  }
}

double SurfaceProjection::getRadius(vtkIdType id_node)
{
  vec3_t x;
  m_FGrid->GetPoint(id_node, x.data());
  projectRestricted(x, id_node);
  vtkIdType id_tri = getProjTriangle(id_node);
  if (id_tri == -1) {
    EG_BUG;
  }
  Triangle &T = m_Triangles[id_tri];
  vec2_t r = T.global3DToLocal2D(x);
  double Ra   = m_Radius[T.idA()];
  double Rb   = m_Radius[T.idB()];
  double Rc   = m_Radius[T.idC()];
  double R    = min(Ra, min(Rb, Rc));
  double Rmax = max(Ra, max(Rb, Rc));
  if (Rmax < 1e90) {
    R = Ra + r[0]*(Rb - Ra) + r[1]*(Rc - Ra);
  }
  return R;
}
