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
#include "laplacesmoother.h"
#include <vtkCellLocator.h>
#include <vtkCharArray.h>
#include <vtkGenericCell.h>
#include "guimainwindow.h"

using namespace GeometryTools;

LaplaceSmoother::LaplaceSmoother() : SurfaceOperation()
{
  DebugLevel = 0;
  setQuickSave(true);
}

void LaplaceSmoother::operate()
{
  QSet<int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  foreach (int bc, bcs) {
    GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(grid);
  }
  UpdatePotentialSnapPoints(false, false);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type" );
  QVector<vtkIdType> smooth_node(grid->GetNumberOfPoints(), false);
  {
    l2g_t nodes = m_Part.getNodes();
    foreach (vtkIdType id_node, nodes) {
      smooth_node[id_node] = true;
    }
  }
  setAllSurfaceCells();
  l2g_t  nodes = m_Part.getNodes();
  g2l_t _nodes = m_Part.getLocalNodes();
  QVector<QVector<int> > n2bc(nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    QSet<int> bcs;
    for (int j = 0; j < m_Part.n2cLSize(i_nodes); ++j) {
      bcs.insert(cell_code->GetValue(m_Part.n2cLG(i_nodes, j)));
    }
    n2bc[i_nodes].resize(bcs.size());
    qCopy(bcs.begin(), bcs.end(), n2bc[i_nodes].begin());
  }
  for (int i_iter = 0; i_iter < m_NumberOfIterations; ++i_iter) {
    foreach (vtkIdType id_node, nodes) {
      if (smooth_node[id_node]) {
        if (node_type->GetValue(id_node) != VTK_FIXED_VERTEX) {
          QVector<vtkIdType> snap_points = getPotentialSnapPoints(id_node);
          if (snap_points.size() > 0) {
            vec3_t x_old;
            vec3_t x_new(0,0,0);
            vec3_t x;
            grid->GetPoint(id_node, x_old.data());
            foreach (vtkIdType id_snap_node, snap_points) {
              grid->GetPoint(id_snap_node, x.data());
              x_new += x;
            }
            x_new *= 1.0/snap_points.size();
            if (n2bc[_nodes[id_node]].size() == 1) {
              int bc = n2bc[_nodes[id_node]][0];
              x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new, id_node);
            } else {
              for (int i_proj_iter = 0; i_proj_iter < 20; ++i_proj_iter) {
                foreach (int bc, n2bc[_nodes[id_node]]) {
                  x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new, id_node);
                }
              }
            }
            grid->GetPoints()->SetPoint(id_node, x_new.data());
          }
        }
      }
    }
  }
}
