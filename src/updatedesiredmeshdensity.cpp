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
#include "updatedesiredmeshdensity.h"
#include "guimainwindow.h"

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity() : SurfaceOperation()
{
  EG_TYPENAME;
  m_MaxEdgeLength = 1e99;
  m_NodesPerQuarterCircle = 0;
  getSet("surface meshing", "minmal number of cells across", 0, m_MinMumCellsAcross);
}


void UpdateDesiredMeshDensity::computeExistingLengths()
{
  QSet<int> all_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(all_bcs);
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
        vec3_t x[N_pts];
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

void UpdateDesiredMeshDensity::operate()
{
  m_ELSManager.read();
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2n   = getPartN2N();
  l2l_t  c2c   = getPartC2C();

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

    // compute node normals
    for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
      normals[i_cells] = GeometryTools::cellNormal(m_Grid, cells[i_cells]);
      normals[i_cells].normalise();
      centres[i_cells] = cellCentre(m_Grid, cells[i_cells]);
    }

    // compute characteristic length according to nodes per quarter circle
    for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
      vec3_t xi = centres[i_cells];
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(cells[i_cells], N_pts, pts);
      for (int j = 0; j < c2c[i_cells].size(); ++j) {
        int j_cells = c2c[i_cells][j];
        if (cell_code->GetValue(cells[i_cells]) == cell_code->GetValue(cells[j_cells])) {
          vec3_t xj = centres[j_cells];
          double cosa  = normals[i_cells]*normals[j_cells];
          double alpha = acos(cosa);
          if (alpha > 0.01*M_PI) {
            vec3_t va = xi - xj;
            vec3_t n1 = normals[i_cells] + normals[j_cells];
            vec3_t n2 = GeometryTools::orthogonalVector(n1);
            n2.normalise();
            va -= (va*n2)*n2;
            double a  = va.abs();
            double R  = 0.5*a/sin(alpha);
            double cl = 0.5*R*M_PI/m_NodesPerQuarterCircle;
            for (int k = 0; k < N_pts; ++k) {
              cl_pre[_nodes[pts[k]]] = min(cl_pre[_nodes[pts[k]]], cl);
            }
          }
        }
      }
    }
  }

  // cells across branches
  computeNormals();
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    QSet<vtkIdType> local_nodes;
    local_nodes.insert(id_node);
    for (int i_level = 0; i_level < 3*m_MinMumCellsAcross; ++i_level) {
      foreach (vtkIdType id_node, local_nodes) {
        for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          local_nodes.insert(m_Part.n2nGG(id_node, i));
        }
      }
    }
    double scal_min = 0.1;
    bool found = false;
    double d = 1e99;
    foreach (vtkIdType id_node1, local_nodes) {
      foreach (vtkIdType id_node2, local_nodes) {
        if (id_node1 != id_node2) {
          vec3_t n1 = m_NodeNormal[id_node1];
          vec3_t n2 = m_NodeNormal[id_node2];
          vec3_t x1, x2;
          m_Grid->GetPoint(id_node1, x1.data());
          m_Grid->GetPoint(id_node2, x2.data());
          vec3_t u = x2 - x1;
          double uabs = u.abs();
          u.normalise();
          double scal1 = n1*n2;
          double scal2 = n1*u;
          if (scal1 < scal_min && scal2 > 0.1) {
            scal_min = scal1;
            found = true;
            d = scal2*uabs;
          }
        }
      }
    }
    if (found) {
      cl_pre[i_nodes] = min(cl_pre[i_nodes], d/m_MinMumCellsAcross);
    }
  }

  // set everything to desired mesh density and find maximal mesh-density
  double cl_min = 1e99;
  int i_nodes_min = -1;
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    double cl = 1e99;
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
    
    if(cl == 0) EG_BUG;
    characteristic_length_desired->SetValue(id_node, cl);
    
    if (cl < cl_min) {
      cl_min = cl;
      i_nodes_min = i_nodes;
    }
  }
  if (i_nodes_min == -1) {
    EG_BUG;
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

  // do a simple averaging step
  /*
  QVector<double> cl_save(m_Grid->GetNumberOfPoints(), 0.0);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    cl_save[id_node] = characteristic_length_desired->GetValue(id_node);
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    double cl_new = 0;
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      cl_new += cl_save[m_Part.n2nGG(id_node, i)];
    }
    if (m_Part.n2nGSize(id_node) > 0) {
      characteristic_length_desired->SetValue(id_node, cl_new/m_Part.n2nGSize(id_node));
    }
  }
  */
}
