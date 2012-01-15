//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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

#include "createhexcore.h"

CreateHexCore::CreateHexCore(vec3_t x1, vec3_t x2, vec3_t xi)
{
  m_X1 = x1;
  m_X2 = x2;
  m_Xi = xi;
}

void CreateHexCore::refineOctree()
{
  m_Octree.resetRefineMarks();
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  do {
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      bool is_surf_node = false;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        if (isSurface(m_Part.n2cGG(id_node, i), m_Grid)) {
          is_surf_node = true;
          break;
        }
      }
      if (is_surf_node) {
        vec3_t x;
        m_Grid->GetPoint(id_node, x.data());
        int i_otcell = m_Octree.findCell(x);
        double h = m_Octree.getDx(i_otcell);
        h = max(h, m_Octree.getDy(i_otcell));
        h = max(h, m_Octree.getDz(i_otcell));
        if (h > cl->GetValue(id_node)) {
          m_Octree.markToRefine(i_otcell);
        }
      }
    }
  } while (m_Octree.refineAll());
}

void CreateHexCore::operate()
{
  m_Octree.setBounds(m_X1, m_X2);
}
