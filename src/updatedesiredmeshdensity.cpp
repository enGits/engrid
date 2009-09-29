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
#include "updatedesiredmeshdensity.h"
#include "guimainwindow.h"

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity() : SurfaceOperation()
{
  EG_TYPENAME;
  m_MaxEdgeLength = 1e99;
  m_NodesPerQuarterCircle = 0;
}


void UpdateDesiredMeshDensity::computeExistingLengths()
{
  QSet<int> all_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(all_bcs);
  QSet<int> fixed_bcs = all_bcs - m_BoundaryCodes;
  QVector<double> edge_length(grid->GetNumberOfPoints(), 0);
  QVector<int> edge_count(grid->GetNumberOfPoints(), 0);
  m_Fixed.fill(false, grid->GetNumberOfPoints());
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      if (fixed_bcs.contains(cell_code->GetValue(id_cell))) {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        vec3_t x[N_pts];
        for (int i = 0; i < N_pts; ++i) {
          grid->GetPoint(pts[i], x[i].data());
          m_Fixed[pts[i]] = true;
        }
        for (int i = 0; i < N_pts; ++i) {
          int j = i + 1;
          if (j >= N_pts) {
            j = 0;
          }
          double L = (x[i] - x[j]).abs();
          edge_length[pts[i]] += L;
          edge_length[pts[j]] += L;
          ++edge_count[pts[i]];
          ++edge_count[pts[j]];
        }
      }
    }
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    if (edge_count[id_node] > 0) {
      double toto = edge_length[id_node]/edge_count[id_node];
      if(toto==0) EG_BUG;
      characteristic_length_desired->SetValue(id_node, toto);
    }
  }
}

void UpdateDesiredMeshDensity::operate()
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2n   = getPartN2N();
  l2l_t  c2c   = getPartC2C();

  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    characteristic_length_specified, grid, "node_specified_density");

  QVector<vec3_t> normals(cells.size(), vec3_t(0,0,0));
  QVector<vec3_t> centres(cells.size(), vec3_t(0,0,0));
  QVector<double> cl_radius(nodes.size(), 1e99);

  computeExistingLengths();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }

  if (m_NodesPerQuarterCircle > 1e-3) {

    // compute node normals
    for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
      normals[i_cells] = GeometryTools::cellNormal(grid, cells[i_cells]);
      normals[i_cells].normalise();
      centres[i_cells] = cellCentre(grid, cells[i_cells]);
    }

    // compute characteristic length according to nodes per quarter circle
    for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
      vec3_t xi = centres[i_cells];
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(cells[i_cells], N_pts, pts);
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
              cl_radius[_nodes[pts[k]]] = min(cl_radius[_nodes[pts[k]]], cl);
            }
          }
        }
      }
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
    cl = min(cl_radius[i_nodes], cl);
    
    double toto = cl;
    if(toto==0) EG_BUG;
    characteristic_length_desired->SetValue(id_node, toto);
    
    qWarning()<<cl<<" < "<<cl_min;
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
        grid->GetPoint(nodes[i_nodes], xi.data());
        for (int j = 0; j < n2n[i_nodes].size(); ++j) {
          int j_nodes = n2n[i_nodes][j];
          double clj = characteristic_length_desired->GetValue(nodes[j_nodes]);
          if (clj > cli && clj > cl_min) {
            vec3_t xj;
            grid->GetPoint(nodes[j_nodes], xj.data());
            ++num_updated;
            double L_new = min(m_MaxEdgeLength, cli * m_GrowthFactor);
            if (!m_Fixed[nodes[j_nodes]]) {
              
              double toto = min(characteristic_length_desired->GetValue(nodes[j_nodes]), L_new);
/*              qDebug()<<"m_MaxEdgeLength="<<m_MaxEdgeLength;
              qDebug()<<"cli="<<cli;
              qDebug()<<"m_GrowthFactor="<<m_GrowthFactor;
              qDebug()<<"characteristic_length_desired->GetValue(nodes[j_nodes])="<<characteristic_length_desired->GetValue(nodes[j_nodes]);
              qDebug()<<"L_new="<<L_new;*/
              if(toto==0) {
                qWarning()<<"m_MaxEdgeLength="<<m_MaxEdgeLength;
                qWarning()<<"cli="<<cli;
                qWarning()<<"m_GrowthFactor="<<m_GrowthFactor;
                qWarning()<<"characteristic_length_desired->GetValue(nodes[j_nodes])="<<characteristic_length_desired->GetValue(nodes[j_nodes]);
                qWarning()<<"L_new="<<L_new;
                EG_BUG;
              }
              characteristic_length_desired->SetValue(nodes[j_nodes], toto);
              
            }
          }
        }
      }
    }
    cl_min *= m_GrowthFactor;
  } while (num_updated > 0);

}
