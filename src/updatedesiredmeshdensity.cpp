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

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity() : SurfaceOperation()
{
  EG_TYPENAME;
  double max_iter;
  getSet("surface meshing", "mesh density maximal number of iterations", 1000.0, max_iter);
  m_MaxIter = double(max_iter);
  getSet("surface meshing", "mesh density convergence limit", 1000, m_ConvLimit);
  m_ConvLimit = 1.0/m_ConvLimit;
  getSet("surface meshing", "cell growth factor", 0, m_GrowthFactor);
}

void UpdateDesiredMeshDensity::operate()
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t nodes  = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  UpdatePotentialSnapPoints(true);
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");
  
  double diff_start = 0;
  double diff = 0;
  bool first = true;
  int iter = 0;
  do {
    diff = 0;
    foreach (vtkIdType id_node, nodes) {
      double md_new = 0;
      double md0 = node_meshdensity_desired->GetValue(id_node);
      int N = n2n[_nodes[id_node]].size();
      foreach (int i_node_neighbour, n2n[_nodes[id_node]]) {
        vtkIdType id_node_neighbour = nodes[i_node_neighbour];
        double Dmd = node_meshdensity_desired->GetValue(id_node_neighbour) - md0;
        /*
        if (Dmd > 0) {
          Dmd = max(0.0, Dmd - m_GrowthFactor*md0);
        } else {
          Dmd = min(0.0, Dmd + m_GrowthFactor*md0);
        }
        */
        md_new += md0 + Dmd;
      }
      md_new /= N;
      VertexMeshDensity nodeVMD = getVMD(id_node);
      int idx = m_VMDvector.indexOf(nodeVMD);
      node_specified_density->SetValue(id_node, idx);
      if ( idx != -1) {
        //md_new = max(md_new, m_VMDvector[idx].density);
        if (md_new > m_VMDvector[idx].density) {
          md_new = m_GrowthFactor*m_VMDvector[idx].density + (1-m_GrowthFactor)*md_new;
        } else {
          md_new = m_VMDvector[idx].density;
        }
      }
      diff = max(diff, fabs(md_new - node_meshdensity_desired->GetValue(id_node)));
      node_meshdensity_desired->SetValue(id_node, md_new);
    }
    ++iter;
    if (first) {
      first = false;
      diff_start = max(1e-30, diff);
    }
    //cout << diff << ',' << diff/diff_start << endl;
  } while ((diff/diff_start > m_ConvLimit) && (iter < m_MaxIter));

  if (iter >= m_MaxIter) {
    cout << "WARNING: Desired mesh density convergence factor has not been reached!" << endl;
  }
}
