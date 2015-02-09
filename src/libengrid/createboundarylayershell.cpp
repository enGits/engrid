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
#include "surfacemeshsmoother.h"
#include "deletecells.h"
#include "stitchholes.h"

CreateBoundaryLayerShell::CreateBoundaryLayerShell()
{
  m_RestGrid      = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_OriginalGrid  = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_PrismaticGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_Success = false;
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
    m_LayerAdjacentNormals[bc] = n0;
    m_LayerAdjacentOrigins[bc] = x0;
  }

  computeBoundaryLayerVectors();
  makeCopy(m_Grid, m_OriginalGrid);
}

QList<vtkIdType> CreateBoundaryLayerShell::correctAdjacentBC(int bc)
{
  cout << "correcting boundary \"" << qPrintable(GuiMainWindow::pointer()->getBC(bc).getName()) << "\"" << endl;

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  double scal_min = -1;
  int    count    = 0;

  SurfaceMeshSmoother smooth;
  smooth.setGrid(m_Grid);
  smooth.useSimpleCentreScheme();

  QList<vtkIdType> cells;
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid)) {
      if (cell_code->GetValue(id_cell) == bc) {
        cells << id_cell;
      }
    }
  }
  CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(bc);
  smooth.setCells(cells);
  smooth.prepareCadInterface(cad);

  QList<vtkIdType> bad_nodes;
  int new_num_bad = m_Grid->GetNumberOfPoints();
  int old_num_bad = 0;
  updateNodeInfo();
  do {
    old_num_bad = new_num_bad;
    scal_min = 1;
    bad_nodes.clear();
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (!m_BoundaryLayerNode[id_node]){
        bool node_has_bc = false;
        bool all_bcs_adjacent = true;
        for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
          int bcn = m_Part.n2bcG(id_node, i);
          if (bcn == bc) {
            node_has_bc = true;
          }
          if (!m_LayerAdjacentBoundaryCodes.contains(bcn)) {
            all_bcs_adjacent = false;
          }
        }
        if (node_has_bc && all_bcs_adjacent) {
          if (m_Part.n2bcGSize(id_node) == 1) {
            vec3_t xs = smooth.smoothNode(id_node);
            xs = cad->snapNode(id_node, xs);
            if (!checkVector(xs)) {
              EG_ERR_RETURN("error while correcting adjacent boundaries");
            }
            m_Grid->GetPoints()->SetPoint(id_node, xs.data());
          } else if (m_Part.n2bcGSize(id_node) == 2 && node_type->GetValue(id_node) != EG_FIXED_VERTEX) {
            int count = 0;
            vec3_t xs(0, 0, 0);
            for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
              vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
              if (m_Part.n2bcGSize(id_neigh) > 1) {
                vec3_t x;
                m_Grid->GetPoint(id_neigh, x.data());
                xs += x;
                ++count;
              }
            }
            if (count < 2) {
              vec3_t x;
              m_Grid->GetPoint(id_node, x.data());
              xs += x;
              ++count;
            }
            xs *= 1.0/count;
            xs = cad->snapToEdge(xs);
            m_Grid->GetPoints()->SetPoint(id_node, xs.data());
          }

          bool node_bad = false;
          for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
            vtkIdType id_cell = m_Part.n2cGG(id_node, i);
            if (isSurface(id_cell, m_Grid)) {
              if (cell_code->GetValue(id_cell) == bc) {
                CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(cell_code->GetValue(id_cell));
                cad->snap(cellCentre(m_Grid, id_cell));
                vec3_t n = cellNormal(m_Grid, id_cell);
                n.normalise();
                double scal = n*cad->getLastNormal();
                scal_min = min(scal_min, scal);
                if (scal < 0.5 && !node_bad) {
                  bad_nodes << id_node;
                  node_bad = true;
                }
              }
            }
          }
        }
      }
    }
    new_num_bad = bad_nodes.size();
    //reduceSurface();
    //cout << "  " << bad_nodes.size() << " node defects" << endl;
    ++count;
  } while (scal_min < 0.5 && count < 1000);
  if (scal_min < 0.5) {
    //writeGrid(grid, "adjacent_bc_failure");
    //return false;
    //EG_ERR_RETURN("failed to correct adjacent surfaces");
  }
  return bad_nodes;
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
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), false);
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    // reset all node types to EG_SIMPLE_VERTEX to trigger a recalculation of the node types
    node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_QUAD) {
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
  QVector<bool> fix(m_Grid->GetNumberOfPoints(), false);
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    // reset all node types to EG_SIMPLE_VERTEX to trigger a recalculation of the node types
    node_type->SetValue(id_node, EG_SIMPLE_VERTEX);
    for (int i = 0; i < part.n2cGSize(id_node); ++i) {
      if (m_Grid->GetCellType(part.n2cGG(id_node, i)) == VTK_QUAD) {
        fix[id_node] = true;
      }
    }
  }
  /*
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
  */
  smooth.fixNodes(fix);
  smooth();
}


void CreateBoundaryLayerShell::operate()
{
  prepare();
  writeBoundaryLayerVectors("blayer");
  createPrismaticGrid();
  m_Success = true;
  m_Part.trackGrid(m_Grid);

  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    QList<vtkIdType> bad_nodes = correctAdjacentBC(bc);
    if (bad_nodes.size() > 0) {
      bool fixable = true;

      QVector<bool> is_bad_cell(m_Grid->GetNumberOfCells(), false);
      foreach (vtkIdType id_node, bad_nodes) {
        if (m_Part.n2bcGSize(id_node) != 1) {
          cout << "node " << id_node << " cannot be fixed." << endl;
          fixable = false;
          break;
        }
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          is_bad_cell[m_Part.n2cGG(id_node, i)] = true;
        }
      }
      if (fixable) {
        QList<vtkIdType> bad_cells;
        EG_FORALL_CELLS (id_cell, m_Grid) {
          if (is_bad_cell[id_cell]) {
            bad_cells << id_cell;
          }
        }

        DeleteCells del;
        del.setGrid(m_Grid);
        del.setCellsToDelete(bad_cells);
        del();
        StitchHoles stitch(bc);
        stitch.setGrid(m_Grid);
        stitch();
      } else {
        m_Success = false;
        cout << "adjacent patch cannot be corrected!" << endl;
      }

    }
  }
  SwapTriangles swap;
  swap.setGrid(m_Grid);
  QSet<int> swap_codes = getAllBoundaryCodes(m_Grid);
  swap_codes -= m_LayerAdjacentBoundaryCodes;
  swap.setBoundaryCodes(swap_codes);
  swap.setVerboseOff();

  for (int iter = 0; iter < 5; ++iter) {
    cout << "correcting adjacent boundaries\n" << "  iteration: " << iter + 1 << endl;
    swap();
    smoothSurface();
    reduceSurface();
    swap();
  }

}

