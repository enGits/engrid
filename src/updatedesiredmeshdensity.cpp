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
  getSet("surface meshing", "cell size increase factor", 1.5, m_StretchFactor);
}

double UpdateDesiredMeshDensity::computeAverage(vtkIdType id_node)
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_density = 0;
  double avg_density = 0;
  EG_VTKDCN (vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  int N = n2n[_nodes[id_node]].size();
  double md0 = node_meshdensity_desired->GetValue(id_node);
  foreach (int i_node_neighbour, n2n[_nodes[id_node]]) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    double Dmd = node_meshdensity_desired->GetValue(id_node_neighbour) - md0;
    if (Dmd > 0) {
      Dmd /= m_StretchFactor;
    } else {
      Dmd *= m_StretchFactor;
    }
    total_density += md0 + Dmd;
  }
  avg_density = total_density / ( double )N;
  return( avg_density );
}

void UpdateDesiredMeshDensity::operate()
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t nodes = getPartNodes();
  
  //UpdateNodeType();
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");
  
  double diff=Convergence_meshdensity+1;
  bool first = true;
  int iter = 0;
  do {
    first = true;
    foreach (vtkIdType node, nodes) {
      VertexMeshDensity nodeVMD = getVMD(node);
      int idx = VMDvector.indexOf(nodeVMD);
      node_specified_density->SetValue(node, idx);
      if ( idx != -1) { //specified
        node_meshdensity_desired->SetValue(node, VMDvector[idx].density);
      } else { //unspecified
        double D = computeAverage(node);
        if(first) {
          diff = fabs(D - node_meshdensity_desired->GetValue(node));
          first = false;
        } else {
          diff = max(fabs(D - node_meshdensity_desired->GetValue(node)), diff);
        }
        node_meshdensity_desired->SetValue(node, D);
      }
    }
    iter++;
  } while (diff > Convergence_meshdensity && !first && iter < MaxiterDensity);// if first=true, it means no new mesh density has been defined (all densities specified)
  if (iter >= MaxiterDensity) {
    cout << "WARNING: Desired mesh density convergence factor has not been reached!" << endl;
  }
}
