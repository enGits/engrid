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

CreateHexCore::CreateHexCore(vec3_t x1, vec3_t x2, vec3_t xi, int num_inital_refinement_levels)
{
  m_X1 = x1;
  m_X2 = x2;
  m_Xi = xi;
  m_NumInitialRefinementLevels = num_inital_refinement_levels;
  m_NumBreakOutLayers = 1;
}

void CreateHexCore::refineOctree()
{
  m_Octree.resetRefineMarks();
  m_Octree.setSmoothTransitionOn();

  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  for (int i = 0; i < m_NumInitialRefinementLevels; ++i) {
    m_Octree.markAllToRefine();
    m_Octree.refineAll();
  }
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

  for (int layer = 0; layer < m_NumBreakOutLayers; ++layer) {
    for (int i = 0; i < m_Octree.getNumCells(); ++i) {
      if (delete_cell[i] && !m_Octree.hasChildren(i)) {
        for (int j = 0; j < 8; ++j) {
          delete_node[m_Octree.getNode(i,j)] = true;
        }
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
  }
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  m_Octree.toVtkGridPolyhedral(otgrid, true);
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
  m_Part.setAllCells();
  deleteOutside(m_Grid);
}

void CreateHexCore::deleteOutside(vtkUnstructuredGrid *grid)
{
  MeshPartition part(grid, true);
  QVector<bool> is_inside(grid->GetNumberOfCells(), false);
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

void CreateHexCore::createBoundaryFaces()
{
  EG_VTKSP(vtkUnstructuredGrid, new_grid);

  // find all polygons which need to be triangulated and collect the new triangles
  m_Part.setAllCells();
  QList<QVector<vtkIdType> > new_triangles;
  QVector<bool> adapt_cell(m_Grid->GetNumberOfCells(), false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vec3_t xc = cellCentre(m_Grid, id_cell);
    if (isVolume(id_cell, m_Grid)) {
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        if (m_Part.c2cGG(id_cell, i) == -1) {
          QVector<vtkIdType> face;
          getFaceOfCell(m_Grid, id_cell, i, face);
          QVector<QVector<vtkIdType> > triangles;
          triangulatePolygon(m_Grid, face, triangles);
          foreach (QVector<vtkIdType> triangle, triangles) {
            vec3_t x1, x2, x3;
            m_Grid->GetPoint(triangle[0], x1.data());
            m_Grid->GetPoint(triangle[1], x2.data());
            m_Grid->GetPoint(triangle[2], x3.data());
            vec3_t xt = (1.0/3.0)*(x1 + x2 + x3);
            vec3_t nt = GeometryTools::triNormal(x1, x2, x3);
            if (nt*(xt - xc) < 0) {
              swap(triangle[0], triangle[1]);
            }
            new_triangles.append(triangle);
          }
          if (face.size() > 3) {
            adapt_cell[id_cell] = true;
          }
        }
      }
    }
  }

  // allocate memory for the new grid
  allocateGrid(new_grid, m_Grid->GetNumberOfCells() + new_triangles.size(), m_Grid->GetNumberOfPoints());

  // copy nodes
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_node, x.data());
  }

  // copy existing cells and update if required
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType id_new_cell;
    if (!adapt_cell[id_cell]) {
      id_new_cell = copyCell(m_Grid, id_cell, new_grid);
    } else {
      QList<QVector<vtkIdType> > faces_of_cell;
      for (int i = 0; i < m_Part.c2cGSize(id_cell); ++i) {
        QVector<vtkIdType> face;
        getFaceOfCell(m_Grid, id_cell, i, face);
        if (m_Part.c2cGG(id_cell, i) == -1) {
          QVector<QVector<vtkIdType> > triangles;
          triangulatePolygon(m_Grid, face, triangles);
          foreach (QVector<vtkIdType> triangle, triangles) {
            faces_of_cell.append(triangle);
          }
        } else {
          faces_of_cell.append(face);
        }
      }
      EG_VTKSP(vtkIdList, stream);
      int stream_size = 1;
      foreach (QVector<vtkIdType> face, faces_of_cell) {
        stream_size += face.size() + 1;
      }
      stream->SetNumberOfIds(stream_size);
      int id = 0;
      stream->SetId(id++, faces_of_cell.size());
      foreach (QVector<vtkIdType> face, faces_of_cell) {
        stream->SetId(id++, face.size());
        foreach (vtkIdType id_node, face) {
          stream->SetId(id++, id_node);
        }
      }
      id_new_cell = new_grid->InsertNextCell(VTK_POLYHEDRON, stream);
    }
    copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
  }

  // determine new boundary code
  EG_VTKDCC(vtkIntArray, cell_code, new_grid, "cell_code");
  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  int bc_new = 1;
  foreach (int bc, bcs) {
    bc_new = max(bc, bc_new);
  }
  ++bc_new;

  // create triangular boundary faces
  EG_VTKDCC(vtkIntArray, cell_orgdir, new_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, new_grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, new_grid, "cell_voldir");
  GuiMainWindow::pointer()->addBC(bc_new, BoundaryCondition("HexCore", "unknown"));
  foreach (QVector<vtkIdType> face, new_triangles) {
    vtkIdType id_new_face;
    if (face.size() == 3) {
      id_new_face = new_grid->InsertNextCell(VTK_TRIANGLE, 3, face.data());
    } else if (face.size() == 4) {
      id_new_face = new_grid->InsertNextCell(VTK_QUAD, 4, face.data());
    } else {
      id_new_face = new_grid->InsertNextCell(VTK_POLYGON, face.size(), face.data());
    }
    cell_code->SetValue(id_new_face, bc_new);
    cell_orgdir->SetValue(id_new_face, 0);
    cell_curdir->SetValue(id_new_face, 0);
    cell_voldir->SetValue(id_new_face, 0);
  }

  makeCopy(new_grid, m_Grid);
}

void CreateHexCore::operate()
{
  m_Octree.setBounds(m_X1, m_X2, 1, 1, 1);
  refineOctree();
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  transferOctreeGrid();
  createBoundaryFaces();
  UpdateCellIndex(m_Grid);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
}
