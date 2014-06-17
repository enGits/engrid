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
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid) && cell_code->GetValue(id_cell) == bc) {
        vec3_t x = cellCentre(m_Grid, id_cell);
        double l = fabs((x - x0)*n0);
        if (l > 0.1*L) {
          BoundaryCondition boundary_condition = GuiMainWindow::pointer()->getBC(bc);
          QString err_msg = "The boundary \"" + boundary_condition.getName() + "\" is not planar.";
          EG_ERR_RETURN(err_msg);
        }
      }
    }
    m_LayerAdjacentNormals[bc] = n0;
    m_LayerAdjacentOrigins[bc] = x0;
  }

  computeBoundaryLayerVectors();
  makeCopy(m_Grid, m_OriginalGrid);
}

void CreateBoundaryLayerShell::correctAdjacentBC(int bc)
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  vec3_t n0 = m_LayerAdjacentNormals[bc];
  vec3_t x0 = m_LayerAdjacentOrigins[bc];
  double scal_min = -1;
  int count = 0;
  while (scal_min < 0.5 && count < 20) {
    scal_min = 1;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_Part.n2bcGSize(id_node) == 1) {
        if (m_Part.n2bcG(id_node, 0) == bc) {
          vec3_t x(0,0,0);
          for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
            vec3_t xn;
            m_Grid->GetPoint(m_Part.n2nGG(id_node, i), xn.data());
            x += xn;
          }
          x *= 1.0/m_Part.n2nGSize(id_node);
          x -= ((x - x0)*n0)*n0;
          m_Grid->GetPoints()->SetPoint(id_node, x.data());
          for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
            vtkIdType id_cell = m_Part.n2cGG(id_node, i);
            if (isSurface(id_cell, m_Grid)) {
              if (cell_code->GetValue(id_cell) == bc) {
                vec3_t n = cellNormal(m_Grid, id_cell);
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
}

void CreateBoundaryLayerShell::createLayerNodes(vtkIdType id_node)
{
  vec3_t x1;
  m_Grid->GetPoint(id_node, x1.data());
  m_PrismaticGrid->GetPoints()->SetPoint(m_ShellNodeMap[id_node], x1.data());
  vec3_t x2 = x1 + m_BoundaryLayerVectors[id_node];
  vec3_t dx = (1.0/m_NumLayers)*(x2 - x1); // temp solution -> CHANGE!!
  vec3_t x = x1;
  for (int i = 0; i <= m_NumLayers; ++i) {
    m_PrismaticGrid->GetPoints()->SetPoint(i*m_ShellPart.getNumberOfNodes() + m_ShellNodeMap[id_node], x.data());
    x += dx; // temp solution -> CHANGE!!
  }
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

  QVector<QList<int> > n2bc(m_PrismaticGrid->GetNumberOfPoints());

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_ShellNodeMap[id_node] != -1) {
      createLayerNodes(id_node);
      for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
        n2bc[m_ShellNodeMap[id_node]].append(m_Part.n2bcG(id_node, i));
        n2bc[m_ShellNodeMap[id_node] + m_ShellPart.getNumberOfNodes()].append(m_Part.n2bcG(id_node, i));
      }
    }
  }

  EG_VTKDCC(vtkIntArray, cell_code_src, m_Grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_code_dst, m_PrismaticGrid, "cell_code");

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
      QVector<vtkIdType> edge(2);
      edge[0] = m_ShellNodeMap[pts[i_pts]];
      edge[1] = m_ShellNodeMap[pts[0]];
      if (i_pts < 2) {
        edge[1] = m_ShellNodeMap[pts[i_pts+1]];
      }
      QSet<int> edge_codes = m_LayerAdjacentBoundaryCodes;
      qcontIntersection(edge_codes, n2bc[edge[0]], edge_codes);
      qcontIntersection(edge_codes, n2bc[edge[1]], edge_codes);
      if (edge_codes.size() == 1) {
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

  /*
  EG_VTKSP(vtkUnstructuredGrid, noquad_grid);
  makeCopy(m_PrismaticGrid, noquad_grid);
  allocateGrid(m_PrismaticGrid, m_PrismaticGrid->GetNumberOfCells() + m_NumLayers*adjacent_edges.size(), m_PrismaticGrid->GetNumberOfPoints());
  */

  /*
  MeshPartition first_prism_part(m_PrimaticGrid, true);
  m_Part.addPartition(first_prism_part);
  */
}

void CreateBoundaryLayerShell::finalise()
{
}

void CreateBoundaryLayerShell::operate()
{
  prepare();
  //writeBoundaryLayerVectors("blayer");
  createPrismaticGrid();

  foreach (int bc, m_LayerAdjacentBoundaryCodes) {
    correctAdjacentBC(bc);
  }


  // ...

  //finalise();
}

