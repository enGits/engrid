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
  m_UseProjection = true;
  m_UseNormalCorrection = false;
}

bool LaplaceSmoother::setNewPosition(vtkIdType id_node, vec3_t x_new)
{
  using namespace GeometryTools;

  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());
  grid->GetPoints()->SetPoint(id_node, x_new.data());
  bool move = true;

  vec3_t n(0,0,0);
  QVector<vec3_t> cell_normals(m_Part.n2cGSize(id_node));
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    cell_normals[i] = GeometryTools::cellNormal(grid, m_Part.n2cGG(id_node, i));
    n += cell_normals[i];
    cell_normals[i].normalise();
  }
  vec3_t x_summit = x_old + n;
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    vec3_t x[3];
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(m_Part.n2cGG(id_node, i), N_pts, pts);
    if (N_pts != 3) {
      EG_BUG;
    }
    for (int j = 0; j < N_pts; ++j) {
      grid->GetPoint(pts[j], x[j].data());
    }
    if (GeometryTools::tetraVol(x[0], x[1], x[2], x_summit, false) <= 0) {
      move = false;
      break;
    }
  }
  if (move) {
    for (int i = 0; i < cell_normals.size(); ++i) {
      for (int j = 0; j < cell_normals.size(); ++j) {
        if (cell_normals[i]*cell_normals[j] < 0.1) {
          move = false;
          break;
        }
      }
    }
  }

  if (!move) {
    grid->GetPoints()->SetPoint(id_node, x_old.data());
  }
  return move;
}

bool LaplaceSmoother::moveNode(vtkIdType id_node, vec3_t &Dx)
{
  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());
  bool moved = false;
  for (int i_relaxation = 0; i_relaxation < 5; ++i_relaxation) {
    if (setNewPosition(id_node, x_old + Dx)) {
      moved = true;
      break;
    }
    Dx *= 0.5;
  }
  return moved;
}


void LaplaceSmoother::operate()
{
  QSet<int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  if (m_UseProjection) {
    foreach (int bc, bcs) {
      GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(grid);
    }
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

  QVector<vec3_t> x_new(nodes.size());

  for (int i_iter = 0; i_iter < m_NumberOfIterations; ++i_iter) {

    QVector<vec3_t> node_normals;
    if (m_UseNormalCorrection) {
      node_normals.fill(vec3_t(0,0,0), grid->GetNumberOfPoints());
      vec3_t n(0,0,0);
      for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          node_normals[id_node] += GeometryTools::cellNormal(grid, m_Part.n2cGG(id_node, i));
        }
        node_normals[id_node].normalise();
      }
    }

    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      vtkIdType id_node = nodes[i_nodes];
      if (smooth_node[id_node]) {
        if (node_type->GetValue(id_node) != VTK_FIXED_VERTEX) {
          QVector<vtkIdType> snap_points = getPotentialSnapPoints(id_node);
          vec3_t n(0,0,0);
          if (snap_points.size() > 0) {
            vec3_t x_old;
            vec3_t x;
            x_new[i_nodes] = vec3_t(0,0,0);
            grid->GetPoint(id_node, x_old.data());
            foreach (vtkIdType id_snap_node, snap_points) {
              grid->GetPoint(id_snap_node, x.data());
              x_new[i_nodes] += x;
              n += node_normals[id_snap_node];
            }
            n.normalise();
            x_new[i_nodes] *= 1.0/snap_points.size();
            if (m_UseNormalCorrection) {
              vec3_t dx = x_new[i_nodes] - x_old;
              dx = (dx*n)*n;
              x_new[i_nodes] -= dx;
            }
            if (m_UseProjection) {
              if (n2bc[_nodes[id_node]].size() == 1) {
                int bc = n2bc[_nodes[id_node]][0];
                x_new[i_nodes] = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new[i_nodes], id_node);
              } else {
                for (int i_proj_iter = 0; i_proj_iter < 20; ++i_proj_iter) {
                  foreach (int bc, n2bc[_nodes[id_node]]) {
                    x_new[i_nodes] = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new[i_nodes], id_node);
                  }
                }
              }
            }
            vec3_t Dx = x_new[i_nodes] - x_old;
            if (moveNode(id_node, Dx)) {
              x_new[i_nodes] = x_old + Dx;
              grid->GetPoints()->SetPoint(id_node, x_old.data());
            } else {
              x_new[i_nodes] = x_old;
            };
          }
        }
      }
    }
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      vtkIdType id_node = nodes[i_nodes];
      if (smooth_node[id_node]) {
        if (node_type->GetValue(id_node) != VTK_FIXED_VERTEX) {
          grid->GetPoints()->SetPoint(id_node, x_new[i_nodes].data());
        }
      }
    }
  }
}
