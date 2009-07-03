//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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

SurfaceProjection::SurfaceProjection()
{
  m_BGrid = vtkUnstructuredGrid::New();
  m_Relax = 0.1;
}

void SurfaceProjection::setBackgroundGrid_initOctree()
{
  double bounds[6];
  m_BGrid->GetBounds(bounds);
  vec3_t x1(bounds[0], bounds[2], bounds[4]);
  vec3_t x2(bounds[1], bounds[3], bounds[5]);
  vec3_t xm = 0.5*(x1 + x2);
  double Dx = 0.5*(x2[0]-x1[0]);
  double Dy = 0.5*(x2[1]-x1[1]);
  double Dz = 0.5*(x2[2]-x1[2]);
  double D = max(Dx, max(Dy, Dz));
  x1 = xm - 2*vec3_t(D,D,D);
  x2 = xm + 2*vec3_t(D,D,D);
  m_OTGrid.setBounds(x1, x2);
  m_OTGrid.setSmoothTransitionOff();
}

void SurfaceProjection::setBackgroundGrid_refineFromNodes()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    foreach (vtkIdType id_node, m_Nodes) {
      vec3_t x;
      m_BGrid->GetPoints()->GetPoint(id_node, x.data());
      int i_otcell = m_OTGrid.findCell(x);
      double Dx = m_OTGrid.getDx(i_otcell);
      double Dy = m_OTGrid.getDy(i_otcell);
      double Dz = m_OTGrid.getDz(i_otcell);
      double D = max(Dx, max(Dy, Dz));
      if (D > m_EdgeLength[id_node]) {
        m_OTGrid.markToRefine(i_otcell);
      }
    }
    num_refine = m_OTGrid.refineAll();
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_refineFromEdges()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    for (int i_nodes = 0; i_nodes < m_Nodes.size(); ++i_nodes) {
      vec3_t x1;
      m_BGrid->GetPoints()->GetPoint(m_Nodes[i_nodes], x1.data());
      for (int i_neigh = 0; i_neigh < m_N2N[i_nodes].size(); ++i_neigh) {
        if (i_nodes < i_neigh) {
          vec3_t x2;
          m_BGrid->GetPoints()->GetPoint(m_Nodes[i_neigh], x2.data());
          for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
            if (!m_OTGrid.hasChildren(i_cells)) {
              double Dx = m_OTGrid.getDx(i_cells);
              double Dy = m_OTGrid.getDy(i_cells);
              double Dz = m_OTGrid.getDz(i_cells);
              double D = max(Dx, max(Dy, Dz));
              for (int i_faces = 0; i_faces < 6; ++i_faces) {
                double k;
                if (m_OTGrid.intersectsFace(i_cells, i_faces, x1, x2, k)) {
                  double L = (1-k)*m_EdgeLength[i_nodes] + k*m_EdgeLength[i_neigh];
                  if (D > L) {
                    m_OTGrid.markToRefine(i_cells);
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_refineFromFaces()
{
  int num_refine;
  do {
    num_refine = 0;
    m_OTGrid.resetRefineMarks();
    for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
      vtkIdType Npts, *pts;
      m_BGrid->GetCellPoints(id_cell, Npts, pts);
      if (Npts == 3) {
        vec3_t a, b, c;
        m_BGrid->GetPoints()->GetPoint(pts[0], a.data());
        m_BGrid->GetPoints()->GetPoint(pts[1], b.data());
        m_BGrid->GetPoints()->GetPoint(pts[2], c.data());
        double La = m_EdgeLength[pts[0]];
        double Lb = m_EdgeLength[pts[1]];
        double Lc = m_EdgeLength[pts[2]];
        for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
          if (!m_OTGrid.hasChildren(i_cells) && !m_OTGrid.markedForRefine(i_cells)) {
            QVector<SortedPair<int> > edges;
            m_OTGrid.getEdges(i_cells, edges);
            foreach (SortedPair<int> edge, edges) {
              vec3_t xi;
              vec3_t ri;
              vec3_t x1 = m_OTGrid.getNodePosition(edge.v1);
              vec3_t x2 = m_OTGrid.getNodePosition(edge.v2);
              if (GeometryTools::intersectEdgeAndTriangle(a, b, c, x1, x2, xi, ri)) {
                double L = La + ri[0]*(Lb-La) + ri[1]*(Lc-La);
                double Dx = m_OTGrid.getDx(i_cells);
                double Dy = m_OTGrid.getDy(i_cells);
                double Dz = m_OTGrid.getDz(i_cells);
                double D = max(Dx, max(Dy, Dz));
                if (D > L) {
                  m_OTGrid.markToRefine(i_cells);
                  break;
                }
              }
            }
          }
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
  } while (num_refine > 0);
}

void SurfaceProjection::setBackgroundGrid_computeLevelSet()
{
  QVector<Triangle> triangles(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType Npts, *pts;
    m_BGrid->GetCellPoints(id_cell, Npts, pts);
    if (Npts == 3) {
      m_BGrid->GetPoints()->GetPoint(pts[0], triangles[id_cell].a.data());
      m_BGrid->GetPoints()->GetPoint(pts[1], triangles[id_cell].b.data());
      m_BGrid->GetPoints()->GetPoint(pts[2], triangles[id_cell].c.data());
      triangles[id_cell].g1 = triangles[id_cell].b - triangles[id_cell].a;
      triangles[id_cell].g2 = triangles[id_cell].c - triangles[id_cell].a;
      triangles[id_cell].g3 = triangles[id_cell].g1.cross(triangles[id_cell].g2);
      triangles[id_cell].g3.normalise();
      triangles[id_cell].G.column(0, triangles[id_cell].g1);
      triangles[id_cell].G.column(1, triangles[id_cell].g2);
      triangles[id_cell].G.column(2, triangles[id_cell].g3);
      triangles[id_cell].GI = triangles[id_cell].G.inverse();
    } else {
      EG_ERR_RETURN("only triangles allowed at the moment");
    }
  }
  m_G.fill(1e99, m_OTGrid.getNumNodes());
  for (int i_nodes = 0; i_nodes < m_OTGrid.getNumNodes(); ++i_nodes) {
    vec3_t xp = m_OTGrid.getNodePosition(i_nodes);
    if ((xp - vec3_t(0,0,0.5)).abs() < 0.05) {
      cout << "here" << endl;
    }
    int pass = 0;
    double total_weight = 0;
    do {
      ++pass;
      foreach (Triangle T, triangles) {
        vec3_t xi;
        vec3_t ri;
        double scal = (xp - T.a)*T.g3;
        double sign = 1;
        if (scal < 0) {
          sign = -1;
        }
        vec3_t x1, x2;
        if (scal > 0) {
          x1 = xp + T.g3;
          x2 = xp - scal*T.g3 - T.g3;
        } else {
          x1 = xp - T.g3;
          x2 = xp - scal*T.g3 + T.g3;
        }
        if (GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri)) {
          pass = 3;
          double G = (xp-T.a)*T.g3;
          if (fabs(G) < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = G;
          }
        } else if (pass == 2) {
        }
      }
    } while (pass < 2);
  }
}

vec3_t SurfaceProjection::calcGradG(vec3_t x)
{
  int cell = m_OTGrid.findCell(x);
  vec3_t DG(0,0,0);
  DG[0] += m_G[m_OTGrid.getNode(cell, 1)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 3)] - m_G[m_OTGrid.getNode(cell, 2)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 5)] - m_G[m_OTGrid.getNode(cell, 4)];
  DG[0] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 6)];
  DG[0] /= m_OTGrid.getDx(cell);
  DG[1] += m_G[m_OTGrid.getNode(cell, 2)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 3)] - m_G[m_OTGrid.getNode(cell, 1)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 6)] - m_G[m_OTGrid.getNode(cell, 4)];
  DG[1] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 5)];
  DG[1] /= m_OTGrid.getDy(cell);
  DG[2] += m_G[m_OTGrid.getNode(cell, 4)] - m_G[m_OTGrid.getNode(cell, 0)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 5)] - m_G[m_OTGrid.getNode(cell, 1)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 6)] - m_G[m_OTGrid.getNode(cell, 2)];
  DG[2] += m_G[m_OTGrid.getNode(cell, 7)] - m_G[m_OTGrid.getNode(cell, 3)];
  DG[2] /= m_OTGrid.getDz(cell);
  DG *= 0.25;
  return DG;
}

double SurfaceProjection::calcG(vec3_t x)
{
  //cout << m_OTGrid.getNumCells() << ',' << x << endl;
  int cell = m_OTGrid.findCell(x);
  vec3_t r = x - m_OTGrid.getNodePosition(m_OTGrid.getNode(cell, 0));
  double kx = r[0]/m_OTGrid.getDx(cell);
  double ky = r[1]/m_OTGrid.getDy(cell);
  double kz = r[2]/m_OTGrid.getDz(cell);
  double g_01 = (1-kx)*m_G[m_OTGrid.getNode(cell, 0)] + kx*m_G[m_OTGrid.getNode(cell, 1)];
  double g_23 = (1-kx)*m_G[m_OTGrid.getNode(cell, 2)] + kx*m_G[m_OTGrid.getNode(cell, 3)];
  double g_45 = (1-kx)*m_G[m_OTGrid.getNode(cell, 4)] + kx*m_G[m_OTGrid.getNode(cell, 5)];
  double g_67 = (1-kx)*m_G[m_OTGrid.getNode(cell, 6)] + kx*m_G[m_OTGrid.getNode(cell, 7)];
  double g_01_23 = (1-ky)*g_01 + ky*g_23;
  double g_45_67 = (1-ky)*g_45 + ky*g_67;
  return (1-kz)*g_01_23 + kz*g_45_67;
}

void SurfaceProjection::writeOctree(QString file_name)
{
  EG_VTKSP(vtkUnstructuredGrid, otg);
  m_OTGrid.toVtkGrid(otg);
  EG_VTKSP(vtkDoubleArray, g);
  g->SetName("g");
  g->SetNumberOfValues(otg->GetNumberOfPoints());
  otg->GetPointData()->AddArray(g);
  for (int i = 0; i < otg->GetNumberOfPoints(); ++i) {
    g->SetValue(i, m_G[i]);
  }
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName((GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtu").toAscii().data());
  vtu->SetDataModeToBinary();
  vtu->SetInput(otg);
  vtu->Write();
  //writeGrid(m_BGrid, "m_BGrid");
}

vec3_t SurfaceProjection::project(vec3_t x)
{
  int count = 0;
  double& _x = x[0];
  double& _y = x[1];
  double& _z = x[2];
  double g = calcG(x);
  while (fabs(g) > 1e-8) {
    if (count > 100) {
      EG_ERR_RETURN("projection failed");
    }
    ++count;
    vec3_t dx = calcGradG(x);
    dx.normalise();
    if (g != 0) {
      x -= m_Relax*g*dx;
    }
    g = calcG(x);
  }
  return x;
}
