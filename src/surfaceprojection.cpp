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
          double L = min(m_EdgeLength[i_nodes], m_EdgeLength[i_neigh]);
          for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
            if (!m_OTGrid.hasChildren(i_cells)) {
              double Dx = m_OTGrid.getDx(i_cells);
              double Dy = m_OTGrid.getDy(i_cells);
              double Dz = m_OTGrid.getDz(i_cells);
              double D = max(Dx, max(Dy, Dz));
              if (D > L) {
                for (int i_faces = 0; i_faces < 6; ++i_faces) {
                  if (m_OTGrid.intersectsFace(i_cells, i_faces, x1, x2)) {
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
        double L = min(m_EdgeLength[pts[0]], min(m_EdgeLength[pts[1]], m_EdgeLength[pts[2]]));
        for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
          if (!m_OTGrid.hasChildren(i_cells) && !m_OTGrid.markedForRefine(i_cells)) {
            double Dx = m_OTGrid.getDx(i_cells);
            double Dy = m_OTGrid.getDy(i_cells);
            double Dz = m_OTGrid.getDz(i_cells);
            double D = max(Dx, max(Dy, Dz));
            if (D > L) {
              QVector<SortedPair<int> > edges;
              m_OTGrid.getEdges(i_cells, edges);
              foreach (SortedPair<int> edge, edges) {
                vec3_t xi;
                vec3_t x1 = m_OTGrid.getNodePosition(edge.v1);
                vec3_t x2 = m_OTGrid.getNodePosition(edge.v2);
                if (GeometryTools::intersectEdgeAndTriangle(a, b, c, x1, x2, xi)) {
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
    foreach (Triangle T, triangles) {
      vec3_t xi;
      vec3_t xp = m_OTGrid.getNodePosition(i_nodes);
      double scal = (xp - T.a)*T.g3;
      double sign = 1;
      if (scal < 0) {
        sign = -1;
      }
      if (GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, xp, xp - 2*scal*T.g3, xi)) {
        double G = (xp-T.a)*T.g3;
        if (fabs(G) < fabs(m_G[i_nodes])) {
          m_G[i_nodes] = G;
        }
      } else {
        double kab = GeometryTools::intersection(T.a, T.b-T.a, xp, T.b-T.a);
        double kac = GeometryTools::intersection(T.a, T.c-T.a, xp, T.c-T.a);
        double kbc = GeometryTools::intersection(T.b, T.c-T.b, xp, T.c-T.b);
        if ((kab >= 0) && (kab <= 1)) {
          xi = T.a + kab*(T.b-T.a);
          double G = (xi-xp).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
        if ((kac >= 0) && (kac <= 1)) {
          xi = T.a + kac*(T.c-T.a);
          double G = (xi-xp).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
        if ((kbc >= 0) && (kbc <= 1)) {
          xi = T.b + kbc*(T.c-T.b);
          double G = (xi-xp).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
        {
          double G = (xp-T.a).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
        {
          double G = (xp-T.b).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
        {
          double G = (xp-T.c).abs();
          if (G < fabs(m_G[i_nodes])) {
            m_G[i_nodes] = sign*G;
          }
        }
      }
    }
  }
}
