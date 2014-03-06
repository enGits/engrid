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
#include "triangularcadinterface.h"
#include "guimainwindow.h"

vtkIdType TriangularCadInterface::m_LastPindex = 0;

TriangularCadInterface::TriangularCadInterface()
{
  m_BGrid = vtkUnstructuredGrid::New();
  this->setGrid(m_BGrid);
  m_CritDistance = 0.1;
  setName("triangulated geometry CAD interface");
}

void TriangularCadInterface::computeSurfaceCurvature()
{
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
      }
    }
  }

  // compute weighted (distance) average of radii
  QVector<double> R_new(m_BGrid->GetNumberOfPoints(), 1e99);
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    vec3_t x1;
    m_BGrid->GetPoint(id_node, x1.data());
    int N = m_BPart.n2nGSize(id_node);
    QVector<double> L(N);
    QVector<double> R(N, -1);
    double Lmax = 0;
    bool average = false;
    for (int i = 0; i < N; ++i) {
      vtkIdType id_neigh = m_BPart.n2nGG(id_node, i);
      vec3_t x2;
      m_BGrid->GetPoint(id_neigh, x2.data());
      L[i] = (x2 - x1).abs();
      if (m_Radius[id_neigh] < 1e90 && L[i] > 0) {
        R[i] = m_Radius[id_neigh];
        Lmax = max(Lmax, L[i]);
        average = true;
      }
    }
    if (average) {
      R_new[id_node] = 0;
      double total_weight = 0;
      for (int i = 0; i < N; ++i) {
        if (R[i] > 0) {
          R_new[id_node] += Lmax/L[i] * R[i];
          total_weight += Lmax/L[i];
          //R_new[id_node] += R[i];
          //total_weight += 1.0;
        }
      }
      R_new[id_node] /= total_weight;
    }
  }
  m_Radius = R_new;
}

void TriangularCadInterface::updateBackgroundGridInfo()
{
  getAllCells(m_Cells, m_BGrid);
  getNodesFromCells(m_Cells, m_Nodes, m_BGrid);
  setBoundaryCodes(GuiMainWindow::pointer()->getAllBoundaryCodes());
  setAllCells();
  readVMD();
  updatePotentialSnapPoints();
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

void TriangularCadInterface::searchNewTriangle(vec3_t xp, vtkIdType &id_tri, vec3_t &x_proj, vec3_t &r_proj, bool &on_triangle)
{
  x_proj = vec3_t(1e99, 1e99, 1e99);
  r_proj = vec3_t(0, 0, 0);
  vec3_t xi, ri;
  double d_min      = 1e99;
  bool   x_proj_set = false;
  on_triangle = false;
  QVector<vtkIdType> candidate_faces;
  m_FaceFinder.getCloseFaces(xp, candidate_faces);
  candidate_faces.clear();
  if (candidate_faces.size() == 0) {
    // backup method -- really inefficient!
    candidate_faces.resize(m_Triangles.size());
    for (vtkIdType id_triangle = 0; id_triangle < m_Triangles.size(); ++id_triangle) {
      candidate_faces[id_triangle] = id_triangle;
    }
  }
  foreach (vtkIdType id_triangle, candidate_faces) {
    Triangle T = m_Triangles[id_triangle];
    double d;
    int side;
    bool intersects = T.snapOntoTriangle(xp, xi, ri, d, side, true);
    if (d >= 1e99) {
      EG_BUG;
    }
    if (d < d_min && (intersects || !on_triangle)) {
      x_proj = xi;
      x_proj_set = true;
      r_proj = ri;
      d_min = d;
      id_tri = id_triangle;
      on_triangle = intersects;
    }
  }
  if (!x_proj_set) {
    EG_BUG;
  }
}

vec3_t TriangularCadInterface::correctCurvature(vtkIdType proj_triangle, vec3_t x)
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

vec3_t TriangularCadInterface::correctCurvature(vec3_t x)
{
  vec3_t x_proj(1e99, 1e99, 1e99);
  vec3_t r_proj(0, 0, 0);
  vtkIdType proj_triangle;
  bool on_triangle = false;
  searchNewTriangle(x, proj_triangle, x_proj, r_proj, on_triangle);
  if (proj_triangle >= m_Triangles.size()) {
    EG_BUG;
  }
  return correctCurvature(proj_triangle, x)  ;
}

vec3_t TriangularCadInterface::snap(vec3_t x, bool correct_curvature)
{
  vec3_t x_proj(1e99, 1e99, 1e99);
  vec3_t r_proj(0, 0, 0);

  bool on_triangle = false;

  vtkIdType proj_triangle;
  searchNewTriangle(x, proj_triangle, x_proj, r_proj, on_triangle);
  if (proj_triangle >= m_Triangles.size()) {
    EG_BUG;
  }
  Triangle T = m_Triangles[proj_triangle];
  vec3_t xi, ri;
  double d;
  int side;
  T.snapOntoTriangle(x, xi, ri, d, side, true);
  x_proj = xi;
  if (x_proj[0] > 1e98) { // should never happen
    EG_BUG;
  }
  r_proj = ri;
  if (!checkVector(x_proj)) {
    x_proj = x;
    m_Failed = true;
  } else {
    vec2_t r = T.global3DToLocal2D(x_proj);
    double Ra   = m_Radius[T.idA()];
    double Rb   = m_Radius[T.idB()];
    double Rc   = m_Radius[T.idC()];
    double R    = min(Ra, min(Rb, Rc));
    double Rmax = max(Ra, max(Rb, Rc));
    if (Rmax < 1e90) {
       R = Ra + r[0]*(Rb - Ra) + r[1]*(Rc - Ra);
    }
    m_LastRadius = R;
    m_LastNormal = T.g3();
    if (correct_curvature) {
      vec3_t x_corr = correctCurvature(proj_triangle, x_proj);
      if (checkVector(x_corr)) {
        x_proj = x_corr;
      }
    }
    m_Failed = false;
  }
  return x_proj;
}

vtkIdType TriangularCadInterface::getProjTriangle(vtkIdType id_node)
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

void TriangularCadInterface::setProjTriangle(vtkIdType id_node, vtkIdType proj_triangle)
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

vec3_t TriangularCadInterface::snapNode(vtkIdType id_node, vec3_t x, bool correct_curvature)
{
  //return snap(x, correct_curvature);

  if (!checkVector(x)) {
    EG_BUG;
  }

  if (id_node < 0) {
    EG_BUG;
  }
  if (id_node >= m_FGrid->GetNumberOfPoints()) {
    EG_BUG;
  }
  vec3_t x_proj(1e99, 1e99, 1e99);
  vec3_t r_proj(0, 0, 0);
  bool x_proj_set = false;

  bool on_triangle = false;

  vtkIdType proj_triangle = -1;
  if (id_node != -1) {
    //proj_triangle = getProjTriangle(id_node);
  }

  if (proj_triangle == -1) {
    searchNewTriangle(x, proj_triangle, x_proj, r_proj, on_triangle);
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
  bool intersects = T.snapOntoTriangle(x, xi, ri, d, side, true);
  if (!intersects || (d > m_CritDistance*T.smallestLength())) {
    searchNewTriangle(x, proj_triangle, x_proj, r_proj, on_triangle);
    T = m_Triangles[proj_triangle];
    T.snapOntoTriangle(x, xi, ri, d, side, true);
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
    x_proj = x;
    m_Failed = true;
  } else {
    vec2_t r = T.global3DToLocal2D(x_proj);
    double Ra   = m_Radius[T.idA()];
    double Rb   = m_Radius[T.idB()];
    double Rc   = m_Radius[T.idC()];
    double R    = min(Ra, min(Rb, Rc));
    double Rmax = max(Ra, max(Rb, Rc));
    if (Rmax < 1e90) {
       R = Ra + r[0]*(Rb - Ra) + r[1]*(Rc - Ra);
    }
    m_LastRadius = R;
    m_LastNormal = T.g3();
    if (correct_curvature) {
      vec3_t x_corr = correctCurvature(proj_triangle, x_proj);
      if (checkVector(x_corr)) {
        x_proj = x_corr;
      }
    }
    m_Failed = false;
  }
  return x_proj;
}

