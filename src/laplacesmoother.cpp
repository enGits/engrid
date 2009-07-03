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
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  QVector<vtkIdType> smooth_node(grid->GetNumberOfPoints(), false);
  {
    l2g_t nodes = m_Part.getNodes();
    foreach (vtkIdType id_node, nodes) {
      smooth_node[id_node] = true;
    }
  }
  setAllSurfaceCells();
  l2g_t nodes = m_Part.getNodes();
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
    cout << "laplace smoother: " << i_iter+1 << "/" << m_NumberOfIterations << endl;
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if ((n2bc[i_nodes].size() < 3) && smooth_node[nodes[i_nodes]]) {
        vec3_t x_old;
        vec3_t x_new(0,0,0);
        vec3_t x;
        int N = 0;
        grid->GetPoint(nodes[i_nodes], x_old.data());
        for (int j = 0; j < m_Part.n2nLSize(i_nodes); ++j) {
          bool use_node = true;
          int j_nodes = m_Part.n2nLL(i_nodes, j);
          foreach (int bc, n2bc[i_nodes]) {
            if (!n2bc[j_nodes].contains(bc)) {
              use_node = false;
              break;
            }
          }
          if (use_node) {
            grid->GetPoint(nodes[j_nodes], x.data());
            x_new += x;
            ++N;
          }
        }
        if (N > 0) {
          x_new *= 1.0/N;
          for (int i_proj_iter = 0; i_proj_iter < 20; ++i_proj_iter) {
            foreach (int bc, n2bc[i_nodes]) {
              x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new);
            }
          }
          grid->GetPoints()->SetPoint(nodes[i_nodes], x_new.data());
        }
      }
    }
  }
}
