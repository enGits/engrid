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
  getSet("surface meshing", "projection octree resolution", 0.5, m_RelEdgeLength);
  getSet("surface meshing", "projection relaxation", 0.5, m_Relax);
  getSet("surface meshing", "projection distance weighting", 1.0, m_DistWeight);
  getSet("surface meshing", "projection distance exponent", 1.0, m_DistExp);
  getSet("surface meshing", "projection direction weighting", 1.0, m_DirWeight);
  getSet("surface meshing", "projection direction exponent", 1.0, m_DirExp);
  getSet("surface meshing", "projection weight offset", 0.001, m_WeightOffset);
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
      if (D > m_RelEdgeLength*m_EdgeLength[id_node]) {
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
                  if (D > m_RelEdgeLength*L) {
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
                if (D > m_RelEdgeLength*L) {
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
  // create m_Triangles
  m_Triangles.resize(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType Npts, *pts;
    m_BGrid->GetCellPoints(id_cell, Npts, pts);
    if (Npts == 3) {
      m_BGrid->GetPoints()->GetPoint(pts[0], m_Triangles[id_cell].a.data());
      m_BGrid->GetPoints()->GetPoint(pts[1], m_Triangles[id_cell].b.data());
      m_BGrid->GetPoints()->GetPoint(pts[2], m_Triangles[id_cell].c.data());
      m_Triangles[id_cell].id_a = pts[0];
      m_Triangles[id_cell].id_b = pts[1];
      m_Triangles[id_cell].id_c = pts[2];
      m_Triangles[id_cell].g1 = m_Triangles[id_cell].b - m_Triangles[id_cell].a;
      m_Triangles[id_cell].g2 = m_Triangles[id_cell].c - m_Triangles[id_cell].a;
      m_Triangles[id_cell].g3 = m_Triangles[id_cell].g1.cross(m_Triangles[id_cell].g2);
      m_Triangles[id_cell].A  = 0.5*m_Triangles[id_cell].g3.abs();
      m_Triangles[id_cell].g3.normalise();
      m_Triangles[id_cell].G.column(0, m_Triangles[id_cell].g1);
      m_Triangles[id_cell].G.column(1, m_Triangles[id_cell].g2);
      m_Triangles[id_cell].G.column(2, m_Triangles[id_cell].g3);
      m_Triangles[id_cell].GI = m_Triangles[id_cell].G.inverse();
    } else {
      EG_ERR_RETURN("only triangles allowed at the moment");
    }
  }

  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0,0,0);
  }
  foreach (Triangle T, m_Triangles) {
    m_NodeNormals[T.id_a] += T.A*T.g3;
    m_NodeNormals[T.id_b] += T.A*T.g3;
    m_NodeNormals[T.id_c] += T.A*T.g3;
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node].normalise();
  }

  // create m_Edges
  QSet<SortedPair<vtkIdType> > raw_edges;
  foreach (Triangle T, m_Triangles) {
    SortedPair<vtkIdType> Pab(T.id_a, T.id_b);
    SortedPair<vtkIdType> Pac(T.id_a, T.id_c);
    SortedPair<vtkIdType> Pbc(T.id_b, T.id_c);
    raw_edges.insert(Pab);
    raw_edges.insert(Pac);
    raw_edges.insert(Pbc);
  }
  m_Edges.resize(raw_edges.size());
  {
    int i_edges = 0;
    foreach (SortedPair<vtkIdType> P, raw_edges) {
      Edge E;
      m_BGrid->GetPoint(P.v1, E.a.data());
      m_BGrid->GetPoint(P.v2, E.b.data());
      E.na = m_NodeNormals[P.v1];
      E.nb = m_NodeNormals[P.v2];
      E.v = E.b - E.a;
      E.La = m_EdgeLength[P.v1];
      E.Lb = m_EdgeLength[P.v2];
      m_Edges[i_edges] = E;
      ++i_edges;
    }
  }

  // initialise G
  m_G.fill(0, m_OTGrid.getNumNodes());

  for (int i_nodes = 0; i_nodes < m_OTGrid.getNumNodes(); ++i_nodes) {
    double weight = 0;
    vec3_t xp = m_OTGrid.getNodePosition(i_nodes);
    foreach (Triangle T, m_Triangles) {
      vec3_t xi(1e99,1e99,1e99);
      vec3_t ri;
      double scal = (xp - T.a)*T.g3;
      vec3_t x1, x2;
      if (scal > 0) {
        x1 = xp + T.g3;
        x2 = xp - scal*T.g3 - T.g3;
      } else {
        x1 = xp - T.g3;
        x2 = xp - scal*T.g3 + T.g3;
      }
      double d = 1e99;
      bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri);
      if (!intersects_face) {
        double kab = GeometryTools::intersection(T.a, T.b - T.a, xp, T.b - T.a);
        double kac = GeometryTools::intersection(T.a, T.c - T.a, xp, T.c - T.a);
        double kbc = GeometryTools::intersection(T.b, T.c - T.b, xp, T.c - T.b);
        double dab = (T.a + kab*(T.b-T.a) - xp).abs();
        double dac = (T.a + kac*(T.c-T.a) - xp).abs();
        double dbc = (T.b + kbc*(T.c-T.b) - xp).abs();
        bool set = false;
        if ((kab >= 0) && (kab <= 1)) {
          if (dab < d) {
            xi = T.a + kab*(T.b-T.a);
            d = dab;
            set = true;
          }
        }
        if ((kac >= 0) && (kac <= 1)) {
          if (dac < d) {
            xi = T.a + kac*(T.c-T.a);
            d = dac;
            set = true;
          }
        }
        if ((kbc >= 0) && (kbc <= 1)) {
          if (dbc < d) {
            xi = T.b + kbc*(T.c-T.b);
            d = dbc;
            set = true;
          }
        }
        double da = (T.a - xp).abs();
        double db = (T.b - xp).abs();
        double dc = (T.c - xp).abs();
        if (da < d) {
          xi = T.a;
          d = da;
          set = true;
        }
        if (db < d) {
          xi = T.b;
          d = db;
        }
        if (dc < d) {
          xi = T.c;
          d = dc;
          set = true;
        }
        if (!set) {
          EG_BUG;
        }
      }
      if (xi[0] > 1e98) {
        EG_BUG;
      }      
      double L = 0.3333*(m_EdgeLength[T.id_a] + m_EdgeLength[T.id_b] + m_EdgeLength[T.id_c]);
      vec3_t dx = xp - xi;
      double g = dx*T.g3;
      if (intersects_face) {
        d = fabs(g);
      }
      double w = m_WeightOffset;
      w += m_DistWeight*pow(L/max(1e-3*L, d), m_DistExp);
      if (dx.abs() < 1e-6*L) {
        w += m_DirWeight;
      } else {
        dx.normalise();
        w += m_DirWeight*pow(fabs(dx*T.g3), m_DirExp);
      }
      //cout << dx << ',' << T.g3 << ',' << dx*T.g3 << endl;
      //m_G[i_nodes] += w*((T.a - xp)*T.g3);
      m_G[i_nodes] += w*g;
      weight += w;
    }
    m_G[i_nodes] /= weight;
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
