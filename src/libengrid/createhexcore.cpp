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
#include "guimainwindow.h"

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
      vtkIdType id_surf = -1;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          id_surf = id_node;
          break;
        } else {
          vtkIdType cell_type = m_Grid->GetCellType(id_cell);
          if (cell_type == VTK_WEDGE) {
            vtkIdType N_pts, *pts;
            m_Grid->GetCellPoints(id_cell, N_pts, pts);
            if      (pts[3] == id_node) id_surf = pts[0];
            else if (pts[4] == id_node) id_surf = pts[1];
            else if (pts[5] == id_node) id_surf = pts[2];
          }
        }
      }
      if (id_surf != -1) {
        vec3_t x;
        m_Grid->GetPoint(id_node, x.data());
        int i_otcell = m_Octree.findCell(x);
        double h = m_Octree.getDx(i_otcell);
        h = max(h, m_Octree.getDy(i_otcell));
        h = max(h, m_Octree.getDz(i_otcell));
        if (h > cl->GetValue(id_surf)) {
          m_Octree.markToRefine(i_otcell);
        }
      }
    }
  } while (m_Octree.refineAll());
}

void CreateHexCore::operate()
{
  m_Octree.setBounds(m_X1, m_X2);
  refineOctree();
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  m_Octree.toVtkGrid(otgrid, false);
  saveGrid(otgrid, GuiMainWindow::pointer()->getCwd() + "/ot");
}
