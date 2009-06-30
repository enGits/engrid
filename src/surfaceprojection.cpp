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
              vec3_t xc = m_OTGrid.getCellCentre(i_cells);
              if (D > L) {
                for (int i_faces = 0; i_faces < 6; ++i_faces) {
                  if (m_OTGrid.intersectsFace(i_cells, i_faces, x1, x2)) {
                    if ((xc[0] < -0.7) && (xc[1] < -0.7)) {
                      m_OTGrid.intersectsFace(i_cells, i_faces, x1, x2);
                      cout << xc << endl;
                    }
                    vec3_t xf = m_OTGrid.getFaceCentre(i_cells, i_faces);
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
    for (int i_cells = 0; i_cells < m_OTGrid.getNumCells(); ++i_cells) {
      if (!m_OTGrid.hasChildren(i_cells)) {
        double Dx = m_OTGrid.getDx(i_cells);
        double Dy = m_OTGrid.getDy(i_cells);
        double Dz = m_OTGrid.getDz(i_cells);
        double D = max(Dx, max(Dy, Dz));
        vec3_t x[8];
        for (int i = 0; i < 8; ++i) {
          x[i] = m_OTGrid.getNodePosition(i_cells, i);
        }
        for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
          vtkIdType Npts, *pts;
          m_BGrid->GetCellPoints(id_cell, Npts, pts);
          if (Npts == 3) {
            double L = min(m_EdgeLength[pts[0]], min(m_EdgeLength[pts[1]], m_EdgeLength[pts[2]]));
            if (D > L) {
              vec3_t xt[3];
              for (int i = 0; i < 3; ++i) {
                m_BGrid->GetPoints()->GetPoint(pts[i], xt[i].data());
              }
              vec3_t xi;
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[0], x[1], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[0], x[2], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[0], x[4], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[1], x[3], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[1], x[5], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[2], x[3], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[2], x[6], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[3], x[7], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[4], x[5], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[4], x[6], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[5], x[7], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
              if (GeometryTools::intersectEdgeAndTriangle(xt[0], xt[1], xt[2], x[6], x[7], xi)) {
                m_OTGrid.markToRefine(i_cells);
                break;
              }
            }
          }
        }
      }
    }
    num_refine = m_OTGrid.refineAll();
  } while (num_refine > 0);
}

