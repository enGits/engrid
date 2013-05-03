//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                     +
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

TriangularCadInterface::TriangularCadInterface()
{
  m_BGrid = vtkUnstructuredGrid::New();
  this->setGrid(m_BGrid);
  m_CritDistance = 0.1;
}

CadInterface::HitType TriangularCadInterface::shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r)
{
  EG_BUG;
  return Miss;
}

/*
CadInterface::PositionType TriangularCadInterface::position(vec3_t x, vec3_t n)
{
  EG_BUG;
  return Outside;
}
*/

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

