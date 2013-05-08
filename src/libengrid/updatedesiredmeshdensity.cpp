//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                      +
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
#include "updatedesiredmeshdensity.h"
#include "guimainwindow.h"
#include "pointfinder.h"

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity() : SurfaceOperation()
{
  EG_TYPENAME;
  m_MaxEdgeLength = 1e99;
  m_NodesPerQuarterCircle = 0;
  m_OnlySurfaceCells = true;
  
  m_GrowthFactor = 0.0;
  m_MinEdgeLength = 0.0;
  m_MinMumCellsAcross = 0;
  m_FeatureResolution2D = 0;
  m_FeatureResolution3D = 0;
  m_FeatureThresholdAngle = deg2rad(45.0);

  getSet("surface meshing", "minimal number of cells across", 0, m_MinMumCellsAcross);
}

double UpdateDesiredMeshDensity::computeSearchDistance(vtkIdType id_face)
{
  vtkIdType N_pts, *pts;
  m_Grid->GetCellPoints(id_face, N_pts, pts);
  QVector<vec3_t> x(N_pts + 1);
  for (int i = 0; i < N_pts; ++i) {
    m_Grid->GetPoint(pts[i], x[i].data());
  }
  x[N_pts] = x[0];
  double L = 0;
  for (int i = 0; i < N_pts; ++i) {
    L = max(L, (x[i] - x[i+1]).abs());
  }
  return L;
}

void UpdateDesiredMeshDensity::computeExistingLengths()
{
  QSet<int> all_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  QSet<int> fixed_bcs = all_bcs - m_BoundaryCodes;
  QVector<double> edge_length(m_Grid->GetNumberOfPoints(), 1e99);
  QVector<int> edge_count(m_Grid->GetNumberOfPoints(), 0);
  m_Fixed.fill(false, m_Grid->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      if (fixed_bcs.contains(cell_code->GetValue(id_cell))) {
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        QVector<vec3_t> x(N_pts);
        for (int i = 0; i < N_pts; ++i) {
          m_Grid->GetPoint(pts[i], x[i].data());
          m_Fixed[pts[i]] = true;
        }
        for (int i = 0; i < N_pts; ++i) {
          int j = i + 1;
          if (j >= N_pts) {
            j = 0;
          }
          double L = (x[i] - x[j]).abs();
          edge_length[pts[i]] = min(edge_length[pts[i]], L);
          edge_length[pts[j]] = min(edge_length[pts[j]], L);
          ++edge_count[pts[i]];
          ++edge_count[pts[j]];
        }
      }
    }
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   m_Grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (edge_count[id_node] > 0) {
      if (edge_length[id_node] > 1e98) {
        EG_BUG;
      }
      characteristic_length_desired->SetValue(id_node, edge_length[id_node]);
    }
  }
}

void UpdateDesiredMeshDensity::computeFeature(const QList<point_t> points, QVector<double> &cl_pre, double res)
{
  int N = 0;
  QVector<vec3_t> pts(points.size());
  for (int i = 0; i < points.size(); ++i) {
    pts[i] = points[i].x;
  }
  PointFinder pfind;
  pfind.setMaxNumPoints(5000);
  pfind.setPoints(pts);
  for (int i = 0; i < points.size(); ++i) {
    double h = 1e99;
    QVector<int> close_points;
    pfind.getClosePoints(points[i].x, close_points, res*points[i].L);
    foreach (int j, close_points) {
      if (i != j) {
        ++N;
        vec3_t x1 = points[i].x;
        vec3_t x2 = points[j].x;
        vec3_t n1 = points[i].n;
        vec3_t n2 = points[j].n;
        vec3_t v = x2 - x1;
        if (n1*n2 < 0) {
          if (n1*v > 0) {
            if (fabs(GeometryTools::angle(n1, (-1)*n2)) <= m_FeatureThresholdAngle) {
              double l = v.abs()/fabs(n1*n2);
              h = min(l/res, h);
            }
          }
        }
      }
    }
    foreach (int i_points, points[i].idx) {
      cl_pre[i_points] = min(h, cl_pre[i_points]);
    }
  }
}

void UpdateDesiredMeshDensity::computeFeature2D(QVector<double> &cl_pre)
{
  if (m_FeatureResolution2D < 1e-3) {
    return;
  }
  QSet<int> bcs = getAllBoundaryCodes(m_Grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (int bc, bcs) {
    QList<point_t> points;
    for (vtkIdType id_face = 0; id_face < m_Grid->GetNumberOfCells(); ++id_face) {
      if (isSurface(id_face, m_Grid)) {
        if (cell_code->GetValue(id_face) == bc) {
          vec3_t xc = cellCentre(m_Grid, id_face);
          vtkIdType num_pts, *pts;
          m_Grid->GetCellPoints(id_face, num_pts, pts);
          QVector<vec3_t> xn(num_pts + 1);
          for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
            m_Grid->GetPoint(pts[i_pts], xn[i_pts].data());
          }
          xn[num_pts] = xn[0];
          for (int i_neigh = 0; i_neigh < m_Part.c2cGSize(id_face); ++i_neigh) {
            vtkIdType id_neigh = m_Part.c2cGG(id_face, i_neigh);
            if (id_neigh != -1) {
              if (cell_code->GetValue(id_neigh) != bc) {
                point_t P;
                P.x = 0.5*(xn[i_neigh] + xn[i_neigh + 1]);
                P.n = xc - xn[i_neigh];
                vec3_t v = xn[i_neigh + 1] - xn[i_neigh];
                v.normalise();
                P.n -= (P.n*v)*v;
                P.n.normalise();
                for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
                  P.idx.append(m_Part.localNode(pts[i_pts]));
                }
                P.L = computeSearchDistance(id_face);
                points.append(P);
              }
            }
          }
        }
      }
    }
    computeFeature(points, cl_pre, m_FeatureResolution2D);
  }
}

void UpdateDesiredMeshDensity::computeFeature3D(QVector<double> &cl_pre)
{
  if (m_FeatureResolution3D < 1e-3) {
    return;
  }
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QList<point_t> points;
  for (vtkIdType id_face = 0; id_face < m_Grid->GetNumberOfCells(); ++id_face) {
    if (isSurface(id_face, m_Grid)) {
      vtkIdType num_pts, *pts;
      m_Grid->GetCellPoints(id_face, num_pts, pts);
      point_t P;
      P.x = cellCentre(m_Grid, id_face);
      P.n = cellNormal(m_Grid, id_face);
      P.n.normalise();
      P.n *= -1;
      for (int i_pts = 0; i_pts < num_pts; ++i_pts) {
        P.idx.append(m_Part.localNode(pts[i_pts]));
      }
      P.L = computeSearchDistance(id_face);
      points.append(P);
    }
  }
  computeFeature(points, cl_pre, m_FeatureResolution3D);
}


void UpdateDesiredMeshDensity::operate()
{
  m_ELSManager.read();
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  
  if (m_OnlySurfaceCells) {
    setAllSurfaceCells();
  } else {
    setAllCells();
  }
  l2g_t  nodes = getPartNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2n   = getPartN2N();

  /*
  QMap<int, double> feature_res;
  QSet<int> bcs = getAllBoundaryCodes(m_Grid);
  foreach (int bc, bcs) {
    feature_res[bc] = computeFeature2D(bc);
    feature_res[bc] = min(feature_res[bc], computeFeature3D(bc));
  }
  */

  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   m_Grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    characteristic_length_specified, m_Grid, "node_specified_density");

  QVector<vec3_t> normals(cells.size(), vec3_t(0,0,0));
  QVector<vec3_t> centres(cells.size(), vec3_t(0,0,0));
  QVector<double> cl_pre(nodes.size(), 1e99);

  computeExistingLengths();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }

  if (m_NodesPerQuarterCircle > 1e-3) {
    QVector<double> R(nodes.size(), 1e99);
    foreach (vtkIdType id_cell, cells) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      int bc = cell_code->GetValue(id_cell);
      for (int i = 0; i < N_pts; ++i) {
        int i_nodes = m_Part.localNode(pts[i]);
        R[i_nodes] = min(R[i_nodes], fabs(GuiMainWindow::pointer()->getCadInterface(bc)->getRadius(pts[i])));
      }
    }
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if (cl_pre[i_nodes] == 0) {
        EG_BUG;
      }
      cl_pre[i_nodes] = max(m_MinEdgeLength, min(cl_pre[i_nodes], 0.5*R[i_nodes]*M_PI/m_NodesPerQuarterCircle));
    }
  }

  // cells across branches
  computeFeature2D(cl_pre);
  computeFeature3D(cl_pre);

  // set everything to desired mesh density and find maximal mesh-density
  double cl_min = 1e99;
  int i_nodes_min = -1;
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    double cl = m_MaxEdgeLength;
    if (m_BoundaryCodes.size() > 0) {
      int idx = characteristic_length_specified->GetValue(id_node);
      if (idx != -1) {
        if (idx >= m_VMDvector.size()) {
          qWarning()<<"idx="<<idx;
          qWarning()<<"m_VMDvector.size()="<<m_VMDvector.size();
          EG_BUG;
        }
        cl = m_VMDvector[idx].density;
      }
    }
    if (m_Fixed[id_node]) {
      cl = characteristic_length_desired->GetValue(id_node);
    }
    cl = min(cl_pre[i_nodes], cl);
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    double cl_src = m_ELSManager.minEdgeLength(x);
    if (cl_src > 0) {
      cl = min(cl, cl_src);
    }
    
    cl = max(m_MinEdgeLength, cl);

    if(cl == 0) {
      EG_BUG;
    }
    characteristic_length_desired->SetValue(id_node, cl);
    
    if (cl < cl_min) {
      for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
        if (m_BoundaryCodes.contains(m_Part.n2bcG(id_node, i))) {
          cl_min = cl;
          i_nodes_min = i_nodes;
          break;
        }
      }
    }
  }
  if (i_nodes_min == -1) {
    EG_ERR_RETURN("There are no edges that need improving.")
  }

  // start from smallest characteristic length and loop as long as nodes are updated
  int num_updated = 0;

  do {
    num_updated = 0;
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      double cli = characteristic_length_desired->GetValue(nodes[i_nodes]);
      if (cli <= cl_min) {
        vec3_t xi;
        m_Grid->GetPoint(nodes[i_nodes], xi.data());
        for (int j = 0; j < n2n[i_nodes].size(); ++j) {
          int j_nodes = n2n[i_nodes][j];
          double clj = characteristic_length_desired->GetValue(nodes[j_nodes]);
          if (clj > cli && clj > cl_min) {
            vec3_t xj;
            m_Grid->GetPoint(nodes[j_nodes], xj.data());
            ++num_updated;
            double L_new = min(m_MaxEdgeLength, cli * m_GrowthFactor);
            if (!m_Fixed[nodes[j_nodes]]) {
              
              double cl_min = min(characteristic_length_desired->GetValue(nodes[j_nodes]), L_new);
              if(cl_min==0) {
                qWarning()<<"m_MaxEdgeLength="<<m_MaxEdgeLength;
                qWarning()<<"cli="<<cli;
                qWarning()<<"m_GrowthFactor="<<m_GrowthFactor;
                qWarning()<<"characteristic_length_desired->GetValue(nodes[j_nodes])="<<characteristic_length_desired->GetValue(nodes[j_nodes]);
                qWarning()<<"L_new="<<L_new;
                EG_BUG;
              }
              characteristic_length_desired->SetValue(nodes[j_nodes], cl_min);
              
            }
          }
        }
      }
    }
    cl_min *= m_GrowthFactor;
  } while (num_updated > 0);

}
