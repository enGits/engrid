// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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

#include "createboundarylayershell.h"

#include "createvolumemesh.h"
#include "swaptriangles.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "meshpartition.h"
#include "deletevolumegrid.h"
#include "laplacesmoother.h"

CreateBoundaryLayerShell::CreateBoundaryLayerShell()
{
  m_RestGrid      = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_OriginalGrid  = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_PrismaticGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
}

void CreateBoundaryLayerShell::prepare()
{
  m_Part.trackGrid(m_Grid);

  DeleteVolumeGrid delete_volume;
  delete_volume.setGrid(m_Grid);
  delete_volume.setAllCells();
  delete_volume();

  readSettings();
  setAllCells();
  getSurfaceCells(m_BoundaryLayerCodes, layer_cells, m_Grid);

  // fill m_LayerAdjacentBoundaryCodes
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (vtkIdType id_cell, layer_cells) {
    for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
      vtkIdType id_neigh = m_Part.c2cGG(id_cell, i);
      int bc = cell_code->GetValue(id_neigh);
      if (!m_BoundaryLayerCodes.contains(bc)) {
        m_LayerAdjacentBoundaryCodes.insert(bc);
      }
    }
  }

  // compute normals and origins of adjacent planes
  m_LayerAdjacentNormals.clear();
  m_LayerAdjacentOrigins.clear();
  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    double L = EG_LARGE_REAL;
    vec3_t n0(0, 0, 0);
    vec3_t x0(0, 0, 0);
    double total_area = 0;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid) && cell_code->GetValue(id_cell) == bc) {
        vec3_t n = cellNormal(m_Grid, id_cell);
        double A = n.abs();
        total_area += A;
        n0 += n;
        x0 += A*cellCentre(m_Grid, id_cell);
        L = min(L, sqrt(4*A/sqrt(3.0)));
      }
    }
    n0.normalise();
    x0 *= 1.0/total_area;
    /*
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid) && cell_code->GetValue(id_cell) == bc) {
        vec3_t x = cellCentre(m_Grid, id_cell);
        double l = fabs((x - x0)*n0);
        if (l > 0.1*L) {
          BoundaryCondition boundary_condition = GuiMainWindow::pointer()->getBC(bc);
          QString err_msg = "The boundary \"" + boundary_condition.getName() + "\" is not planar.\n";
          QString L_txt, l_txt;
          L_txt.setNum(L);
          l_txt.setNum(l);
          err_msg += "L = " + L_txt + " ,  l = " + l_txt;
          EG_ERR_RETURN(err_msg);
        }
      }
    }*/
    m_LayerAdjacentNormals[bc] = n0;
    m_LayerAdjacentOrigins[bc] = x0;
  }

  computeBoundaryLayerVectors();
  makeCopy(m_Grid, m_OriginalGrid);
}

void CreateBoundaryLayerShell::correctAdjacentBC(int bc, vtkUnstructuredGrid *grid)
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  MeshPartition part(grid, true);
  vec3_t n0 = m_LayerAdjacentNormals[bc];
  vec3_t x0 = m_LayerAdjacentOrigins[bc];
  double scal_min = -1;
  int count = 0;
  while (scal_min < 0.5 && count < 1000) {
    scal_min = 1;
    for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
      if (part.n2bcGSize(id_node) == 1) {
        if (part.n2bcG(id_node, 0) == bc) {
          vec3_t x(0,0,0);
          int count = 0;
          double A_total = 0;
          for (int i = 0; i < part.n2cGSize(id_node); ++i) {
            vtkIdType id_cell = part.n2cGG(id_node, i);
            if (isSurface(id_cell, grid)) {
              if (cell_code->GetValue(id_cell) == bc) {
                vec3_t n = cellNormal(grid, id_cell);
                double A = n.abs();
                n.normalise();
                if (n*n0 > 0.5) {
                  x += A*cellCentre(grid, id_cell);
                  A_total += A;
                  ++count;
                }
              }
            }
          }
          if (count == 0) {
            grid->GetPoint(id_node, x.data());
          } else {
            x *= 1.0/A_total;
          }

          x -= ((x - x0)*n0)*n0;
          grid->GetPoints()->SetPoint(id_node, x.data());
          for (int i = 0; i < part.n2cGSize(id_node); ++i) {
            vtkIdType id_cell = part.n2cGG(id_node, i);
            if (isSurface(id_cell, grid)) {
              if (cell_code->GetValue(id_cell) == bc) {
                vec3_t n = cellNormal(grid, id_cell);
                n.normalise();
                scal_min = min(scal_min, n*n0);
              }
            }
          }
        }
      }
    }
    ++count;
  }
  if (scal_min < 0.5) {
    writeGrid(grid, "adjacent_bc_failure");
    EG_ERR_RETURN("failed to correct adjacent surfaces");
  }
}

void CreateBoundaryLayerShell::createLayerNodes(vtkIdType id_node)
{
  vec3_t x1;
  m_Grid->GetPoint(id_node, x1.data());
  m_PrismaticGrid->GetPoints()->SetPoint(m_ShellNodeMap[id_node], x1.data());
  vec3_t x2 = x1 + m_BoundaryLayerVectors[id_node];

  double H  = m_BoundaryLayerVectors[id_node].abs();
  double h  = H*(1.0 - m_StretchingRatio)/(1.0 - pow(m_StretchingRatio, m_NumLayers));
  vec3_t dx = (1.0/H)*m_BoundaryLayerVectors[id_node];
  vec3_t x  = x1;
  m_PrismaticGrid->GetPoints()->SetPoint(m_ShellNodeMap[id_node], x1.data());
  for (int i = 1; i < m_NumLayers; ++i) {
    x += h*dx;
    h *= m_StretchingRatio;
    m_PrismaticGrid->GetPoints()->SetPoint(i*m_ShellPart.getNumberOfNodes() + m_ShellNodeMap[id_node], x.data());
  }
  m_PrismaticGrid->GetPoints()->SetPoint(m_NumLayers*m_ShellPart.getNumberOfNodes() + m_ShellNodeMap[id_node], x2.data());
  m_Grid->GetPoints()->SetPoint(id_node, x2.data());
}

void CreateBoundaryLayerShell::createPrismaticGrid()
{
  QVector<vtkIdType> original_triangles, shell_triangles;
  getSurfaceCells(m_BoundaryLayerCodes, original_triangles, m_OriginalGrid);
  getSurfaceCells(m_BoundaryLayerCodes, shell_triangles, m_Grid);
  {
    MeshPartition part(m_Grid);
    part.setCells(shell_triangles);
    allocateGrid(m_PrismaticGrid, (m_NumLayers + 1)*part.getNumberOfCells(), (m_NumLayers + 1)*part.getNumberOfNodes());
  }

  m_ShellNodeMap.fill(-1, m_Grid->GetNumberOfPoints());
  m_ShellPart.setGrid(m_Grid);
  m_ShellPart.setCells(shell_triangles);
  for (int i = 0; i < m_ShellPart.getNumberOfNodes(); ++i) {
    m_ShellNodeMap[m_ShellPart.globalNode(i)] = i;
  }

  QVector<QSet<int> > n2bc(m_PrismaticGrid->GetNumberOfPoints());

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_ShellNodeMap[id_node] != -1) {
      createLayerNodes(id_node);
      for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
        n2bc[m_ShellNodeMap[id_node]].insert(m_Part.n2bcG(id_node, i));
        //n2bc[m_ShellNodeMap[id_node] + m_ShellPart.getNumberOfNodes()].insert(m_Part.n2bcG(id_node, i));
      }
    }
  }

  QList<QVector<vtkIdType> > adjacent_edges;

  // create prismatic cells and prepare adjacent quad faces
  //
  foreach (vtkIdType id_cell, shell_triangles) {
    vtkIdType num_pts, *pts;
    m_Grid->GetCellPoints(id_cell, num_pts, pts);
    vtkIdType tri_pts[3], pri_pts[6];
    for (int i_pts = 0; i_pts < 3; ++i_pts) {
      if (m_ShellNodeMap[pts[i_pts]] < 0) {
        EG_BUG;
      }
      if (m_ShellNodeMap[pts[i_pts]] >= m_ShellPart.getNumberOfNodes()) {
        EG_BUG;
      }
      QVector<vtkIdType> edge(4);
      edge[1] = m_ShellNodeMap[pts[i_pts]];
      edge[2] = m_ShellNodeMap[pts[0]];
      if (i_pts < 2) {
        edge[2] = m_ShellNodeMap[pts[i_pts+1]];
      }
      QSet<int> edge_codes = m_LayerAdjacentBoundaryCodes;
      edge_codes.intersect(n2bc[edge[1]]);
      edge_codes.intersect(n2bc[edge[2]]);
      if (edge_codes.size() == 1) {
        edge[0] = *edge_codes.begin();
        edge[3] = id_cell;
        adjacent_edges.append(edge);
      }
      tri_pts[i_pts] = m_ShellNodeMap[pts[i_pts]];
    }
    vtkIdType id_tri = m_PrismaticGrid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
    copyCellData(m_Grid, id_cell, m_PrismaticGrid, id_tri);
    for (int i_layer = 0; i_layer < m_NumLayers; ++i_layer) {
      for (int i_pts = 0; i_pts < 3; ++i_pts) {
        pri_pts[i_pts]     = m_ShellNodeMap[pts[i_pts]] + i_layer*m_ShellPart.getNumberOfNodes();
        pri_pts[i_pts + 3] = m_ShellNodeMap[pts[i_pts]] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
      }
      vtkIdType id_pri = m_PrismaticGrid->InsertNextCell(VTK_WEDGE, 6, pri_pts);
    }
  }

  // create quads on adjacent boundary faces
  //
  EG_VTKSP(vtkUnstructuredGrid, noquad_grid);
  makeCopy(m_PrismaticGrid, noquad_grid);
  allocateGrid(m_PrismaticGrid, m_PrismaticGrid->GetNumberOfCells() + m_NumLayers*adjacent_edges.size(), m_PrismaticGrid->GetNumberOfPoints());
  makeCopyNoAlloc(noquad_grid, m_PrismaticGrid);

  EG_VTKDCC(vtkIntArray, cell_code, m_PrismaticGrid, "cell_code");

  foreach (QVector<vtkIdType> edge, adjacent_edges) {
    vtkIdType qua_pts[4];
    for (int i_layer = 0; i_layer < m_NumLayers; ++i_layer) {
      qua_pts[0] = edge[2] + i_layer*m_ShellPart.getNumberOfNodes();
      qua_pts[1] = edge[1] + i_layer*m_ShellPart.getNumberOfNodes();
      qua_pts[2] = edge[1] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
      qua_pts[3] = edge[2] + (i_layer + 1)*m_ShellPart.getNumberOfNodes();
      vtkIdType id_qua = m_PrismaticGrid->InsertNextCell(VTK_QUAD, 4, qua_pts);
      copyCellData(m_Grid, edge[3], m_PrismaticGrid, id_qua);
      cell_code->SetValue(id_qua, edge[0]);
    }
  }
}

void CreateBoundaryLayerShell::reduceSurface()
{
  RemovePoints remove_points;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  remove_points.setMeshPartition(part);
  remove_points.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);
  remove_points.setUpdatePSPOn();
  remove_points.setThreshold(3);
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), true);

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = false;
      }
    }
  }
  for (int layer = 0; layer < 3; ++layer) {
    QVector<bool> tmp = fix;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!tmp[id_node]) {
        for (int i = 0; i < part.n2nGSize(id_node); ++i) {
          fix[part.n2nGG(id_node, i)] = false;
        }
      }
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = true;
      }
    }
  }

  remove_points.fixNodes(fix);
  remove_points();
}

void CreateBoundaryLayerShell::smoothSurface()
{
  LaplaceSmoother smooth;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  smooth.setMeshPartition(part);
  smooth.setNumberOfIterations(2);
  smooth.setBoundaryCodes(m_LayerAdjacentBoundaryCodes);

  QVector<bool> fix(m_Grid->GetNumberOfPoints(), true);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_WEDGE) {
        fix[id_node] = false;
      }
    }
  }
  for (int layer = 0; layer < 3; ++layer) {
    QVector<bool> tmp = fix;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!tmp[id_node]) {
        for (int i = 0; i < part.n2nGSize(id_node); ++i) {
          fix[part.n2nGG(id_node, i)] = false;
        }
      }
    }
  }

  smooth.fixNodes(fix);
  smooth();
}


void CreateBoundaryLayerShell::operate()
{
  prepare();
  writeBoundaryLayerVectors("blayer");
  createPrismaticGrid();
  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    correctAdjacentBC(bc, m_Grid);
  }
  SwapTriangles swap;
  swap.setGrid(m_Grid);
  QSet<int> swap_codes = getAllBoundaryCodes(m_Grid);
  swap_codes -= m_LayerAdjacentBoundaryCodes;
  swap.setBoundaryCodes(swap_codes);
  swap.setVerboseOff();

  for (int iter = 0; iter < 1; ++iter) {
    swap();
    smoothSurface();
    reduceSurface();
    swap();
  }

}

