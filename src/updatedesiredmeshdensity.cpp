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
  getSet("surface meshing", "cell growth factor", 0, m_GrowthFactor);
}

void UpdateDesiredMeshDensity::operate()
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  setAllSurfaceCells();
  l2g_t nodes = getPartNodes();
  l2l_t n2n   = getPartN2N();

  EG_VTKDCN(vtkDoubleArray, md_desired,   grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    md_specified, grid, "node_specified_density");

  QVector<bool> md_set(nodes.size(), false);
  QVector<bool> md_preset(nodes.size(), false);

  // set everything to desired mesh density and find maximal mesh-density
  double md_max = 0;
  int i_nodes_max = -1;
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    double md = 1e-10;
    int idx = md_specified->GetValue(id_node);
    if (idx != -1) {
      if (idx >= m_VMDvector.size()) {
        EG_BUG;
      }
      md = m_VMDvector[idx].density;
    }
    md_desired->SetValue(id_node, md);
    if (md > md_max) {
      md_max = md;
      i_nodes_max = i_nodes;
    }
  }
  if (i_nodes_max == -1) {
    EG_BUG;
  }
  md_set[i_nodes_max] = true;

  // start from highest mesh density and loop as long as nodes are updated
  int num_updated = 0;

  do {
    num_updated = 0;
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if (md_set[i_nodes]) {
        vec3_t xi;
        grid->GetPoint(nodes[i_nodes], xi.data());
        for (int j = 0; j < n2n[i_nodes].size(); ++j) {
          int j_nodes = n2n[i_nodes][j];
          if (!md_set[j_nodes]) {
            vec3_t xj;
            grid->GetPoint(nodes[j_nodes], xj.data());
            if (!md_preset[j_nodes]) {
              md_preset[j_nodes] = true;
              ++num_updated;
            }
            double L_new = 1.0/md_desired->GetValue(nodes[i_nodes]) + (m_GrowthFactor - 1)*(xi-xj).abs();
            md_desired->SetValue(nodes[j_nodes], max(md_desired->GetValue(nodes[j_nodes]), 1.0/L_new));
          }
        }
      }
    }
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if (md_preset[i_nodes]) {
        md_set[i_nodes] = true;
      }
    }
  } while (num_updated > 0);

}
