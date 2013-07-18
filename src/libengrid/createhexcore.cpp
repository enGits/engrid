// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "pointfinder.h"

CreateHexCore::CreateHexCore(vec3_t x1, vec3_t x2, vec3_t xi, int num_i, int num_j, int num_k)
{
  vec3_t xm = 0.5*(x1 + x2);
  vec3_t Dx = x2 - x1;
  m_X1 = x1;
  m_X2 = x2;
  m_Xi = xi;
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
}

void CreateHexCore::refineOctree()
{
  m_Octree.resetRefineMarks();
  m_Octree.setSmoothTransitionOn();

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

void CreateHexCore::transferOctreeGrid()
{
  cout << m_Grid->GetNumberOfPoints() << endl;
  QVector<bool> delete_node(m_Octree.getNumNodes(), false);
  QVector<bool> delete_cell(m_Octree.getNumCells(), false);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    int i_cell = m_Octree.findCell(x);
    if (m_Octree.hasChildren(i_cell)) {
      EG_BUG;
    }
    delete_cell[i_cell] = true;
    for (int i = 0; i < 8; ++i) {
      delete_node[m_Octree.getNode(i_cell, i)] = true;
    }
  }
  for (int i = 0; i < m_Octree.getNumCells(); ++i) {
    if (!m_Octree.hasChildren(i)) {
      for (int j = 0; j < 8; ++j) {
        if (delete_node[m_Octree.getNode(i,j)]) {
          delete_cell[i] = true;
          break;
        }
      }
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  m_Octree.toVtkGridPolyhedral(otgrid, true);
  //writeGrid(otgrid, "debug_01");
  MeshPartition add_part(otgrid);
  QList<vtkIdType> add_cells;
  for (vtkIdType id_cell = 0; id_cell < otgrid->GetNumberOfCells(); ++id_cell) {
    int i_cell = m_Octree.findCell(cellCentre(otgrid, id_cell));
    if (!delete_cell[i_cell]) {
      add_cells.append(id_cell);
    }
  }
  add_part.setCells(add_cells);
  m_Part.addPartition(add_part);  
  //writeGrid(m_Grid, "debug_02");
  m_Part.setAllCells();
  deleteOutside(m_Grid);
}

void CreateHexCore::deleteOutside(vtkUnstructuredGrid *grid)
{
  MeshPartition part(grid, true);
  QVector<bool> is_inside(grid->GetNumberOfCells(), false);
  is_inside.fill(true); // DEBUG
  vtkIdType id_start = -1;
  double dmin = 1e99;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      vec3_t x = cellCentre(grid, id_cell);
      double d = (x - m_Xi).abs();
      if (d < dmin) {
        dmin = d;
        id_start = id_cell;
      }
    }
  }
  if (id_start == -1) {
    EG_BUG;
  }
  is_inside[id_start] = true;
  bool added = true;
  while (added) {
    added = false;
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      if (is_inside[id_cell]) {
        for (int j = 0; j < part.c2cGSize(id_cell); ++j) {
          vtkIdType id_neigh = part.c2cGG(id_cell, j);
          if (id_neigh >= 0) {
            if (!is_inside[id_neigh]) {
              is_inside[id_neigh] = true;
              added = true;
            }
          }
        }
      }
    }
  }
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      is_inside[id_cell] = true;
    }
    if (is_inside[id_cell]) {
      ++N;
    }
  }
  QVector<vtkIdType> cls(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (is_inside[id_cell]) {
      cls[N] = id_cell;
      ++N;
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, return_grid);
  makeCopy(grid, return_grid, cls);
  makeCopy(return_grid, grid);
}

void CreateHexCore::operate()
{
  m_Octree.setBounds(m_X1, m_X2, m_NumI, m_NumJ, m_NumK);
  refineOctree();
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  transferOctreeGrid();
}
