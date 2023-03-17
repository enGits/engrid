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
// 
#include "tetgenoperation.h"
#include "engrid.h"
#include "guimainwindow.h"
#include "optimisenormalvector.h"

#include <QDir>

TetGenOperation::TetGenOperation()
{
  // The following code will be used once the option to start TetGen in a
  // separate process has been implemented

  /*
  QDir cwd(GuiMainWindow::pointer()->getCwd());
  m_TetGenPath = GuiMainWindow::pointer()->getCwd() + "tetgen";
  QDir tetgen_dir(m_TetGenPath);
  if (!tetgen_dir.exists()) {
    cwd.mkdir("tetgen");
  }
  */
}

void TetGenOperation::copyToTetGen(tetgenio &tgio, vtkUnstructuredGrid *alt_grid)
{
  vtkUnstructuredGrid *grid = alt_grid;
  if (!grid) {
    grid = m_Grid;
  }
  MeshPartition part(grid, true);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, grid, "cell_voldir");
  tgio.initialize();
  tgio.numberofpoints = grid->GetNumberOfPoints();
  tgio.pointlist = new REAL [3*tgio.numberofpoints];
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      tgio.pointlist[3*id_node + i] = x[i];
    }
  }
  QList<vtkIdType> triangles, tetrahedra;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if      (grid->GetCellType(id_cell) == VTK_TRIANGLE) triangles.append(id_cell);
    else if (grid->GetCellType(id_cell) == VTK_TETRA)    tetrahedra.append(id_cell);
    else {
      EG_BUG;
    }
  }
  tgio.numberoffacets = triangles.size();
  tgio.facetmarkerlist = new int [triangles.size()];
  tgio.facetlist = new tetgenio::facet [triangles.size()];
  {
    int i = 0;
    foreach (vtkIdType id_cell, triangles) {
      m_OrgDir = cell_orgdir->GetValue(id_cell);
      m_CurDir = cell_curdir->GetValue(id_cell);
      m_VolDir = cell_voldir->GetValue(id_cell);
      EG_GET_CELL(id_cell, grid);
      tgio.facetlist[i].numberofpolygons = 1;
      tgio.facetlist[i].polygonlist = new tetgenio::polygon[1];
      tgio.facetlist[i].numberofholes = 0;
      tgio.facetlist[i].holelist = NULL;
      tetgenio::polygon &poly = tgio.facetlist[i].polygonlist[0];
      poly.numberofvertices = 3;
      poly.vertexlist = new int[3];
      for (int j = 0; j < num_pts; ++j) {
        poly.vertexlist[j] = pts[j];
      }
      tgio.facetmarkerlist[i] = cell_code->GetValue(id_cell);
      ++i;
    }
  }
  {
    tgio.numberoftetrahedra = tetrahedra.size();
    tgio.tetrahedronlist = new int [tetrahedra.size()*4];
    int i = 0;
    foreach (vtkIdType id_cell, tetrahedra) {
      EG_GET_CELL(id_cell, grid);
      for (int j = 0; j < num_pts; ++j) {
        tgio.tetrahedronlist[i*num_pts + j] = pts[j];
      }
      ++i;
    }

    QVector<bool> provide_mesh_resolution(grid->GetNumberOfPoints(), true);
    for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
      provide_mesh_resolution[id_node] = false;
    }
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, grid)) {
        EG_GET_CELL(id_cell, grid);
        for (int i = 0; i < num_pts; ++i) {
          provide_mesh_resolution[pts[i]] = true;
        }
      }
    }
    for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
      for (int i = 0; i < part.n2nGSize(id_node); ++i) {
        provide_mesh_resolution[part.n2nGG(id_node,i)] = true;
      }
    }

    tgio.numberofpointmtrs = 1;
    tgio.pointmtrlist = new REAL [tgio.numberofpoints];
    for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
      tgio.pointmtrlist[id_node] = 0;
      if (provide_mesh_resolution[id_node]) {
        vec3_t x;
        grid->GetPoint(id_node, x.data());        
        /// @todo look into this
        double l = min(m_MaximalEdgeLength, 0.66*m_ELSManager.minEdgeLength(x)); // apply empirical scaling to match surface resolution
        tgio.pointmtrlist[id_node] = l;
      }
    }
  }
}

void TetGenOperation::copyFromTetGen(tetgenio &tgio)
{
  // build search tree to correct surface orientation
  QVector<vtkIdType> tri_faces;
  getAllCellsOfType(VTK_TRIANGLE, tri_faces, m_Grid);
  QVector<Triangle> triangles(tri_faces.size());
  TriangleTree search_tree;
  {
    int i = 0;
    foreach (vtkIdType id_cell, tri_faces) {
      EG_GET_CELL(id_cell, m_Grid);
      vec3_t a, b, c;
      m_Grid->GetPoint(pts[0], a.data());
      m_Grid->GetPoint(pts[1], b.data());
      m_Grid->GetPoint(pts[2], c.data());
      triangles[i] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      ++i;
    }
  }
  search_tree.rebuild(triangles.begin(), triangles.end());
  search_tree.accelerate_distance_queries();

  EG_VTKSP(vtkUnstructuredGrid, full_grid);
  allocateGrid(full_grid, tgio.numberoftetrahedra + tgio.numberoftrifaces, tgio.numberofpoints);
  EG_VTKDCC(vtkIntArray, cell_code, full_grid, "cell_code");
  EG_VTKDCC(vtkIntArray, cell_orgdir, full_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, full_grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, full_grid, "cell_voldir");
  for (vtkIdType id_node = 0; id_node < tgio.numberofpoints; ++id_node) {
    vec3_t x;
    for (int i = 0; i < 3; ++i) {
      x[i] = tgio.pointlist[3*id_node + i];
    }
    full_grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (int i = 0; i < tgio.numberoftetrahedra; ++i) {
    vtkIdType num_pts = 4, pts[4];
    for (int j = 0; j < num_pts; ++j) {
      pts[j] = tgio.tetrahedronlist[i*num_pts + j];
    }
    vtkIdType id_cell = full_grid->InsertNextCell(VTK_TETRA, num_pts, pts);
    cell_code->SetValue(id_cell, 0);
  }
  for (int i = 0; i < tgio.numberoftrifaces; ++i) {
    vtkIdType num_pts = 3, pts[3];
    for (int j = 0; j < num_pts; ++j) {
      //pts[2-j] = tgio.trifacelist[i*num_pts + j];
      pts[j] = tgio.trifacelist[i*num_pts + j];
    }

    vec3_t x1, x2, x3;
    full_grid->GetPoint(pts[0], x1.data());
    full_grid->GetPoint(pts[1], x2.data());
    full_grid->GetPoint(pts[2], x3.data());
    vec3_t xm = (1.0/3.0)*(x1 + x2 + x3);
    vec3_t n1 = GeometryTools::triNormal(x1, x2, x3);
    TriangleTree::Point_and_primitive_id cp = search_tree.closest_point_and_primitive(Point(xm[0], xm[1], xm[2]));

    Point  pa = cp.second->vertex(0);
    vec3_t xa(pa[0], pa[1], pa[2]);
    Point  pb = cp.second->vertex(1);
    vec3_t xb(pb[0], pb[1], pb[2]);
    Point  pc = cp.second->vertex(2);
    vec3_t xc(pc[0], pc[1], pc[2]);

    vec3_t n2 = GeometryTools::triNormal(xa, xb, xc);
    if (n1*n2 < 0) {
      swap(pts[0], pts[1]);
    }

    vtkIdType id_cell = full_grid->InsertNextCell(VTK_TRIANGLE, num_pts, pts);
    cell_code->SetValue(id_cell, tgio.trifacemarkerlist[i]);
    cell_orgdir->SetValue(id_cell, m_OrgDir);
    cell_curdir->SetValue(id_cell, m_CurDir);
    cell_voldir->SetValue(id_cell, m_VolDir);
  }

  QVector<int> keep_cell(full_grid->GetNumberOfCells(), 0);
  MeshPartition full_part(full_grid, true);
  QList<vtkIdType> seed_cells;
  EG_FORALL_CELLS(id_cell, full_grid) {
    EG_GET_CELL(id_cell, full_grid);
    if (type_cell == VTK_TRIANGLE) {
      QVector<QSet<vtkIdType> > cells(3);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < full_part.n2cGSize(pts[i]); ++j) {
          cells[i].insert(full_part.n2cGG(pts[i], j));
        }
      }
      cells[0].intersect(cells[1]);
      cells[0].intersect(cells[2]);
      vec3_t x_face = cellCentre(full_grid, id_cell);
      vec3_t n_face = cellNormal(full_grid, id_cell);
      foreach (vtkIdType id_neigh, cells[0]) {
        if (full_grid->GetCellType(id_neigh) == VTK_TETRA) {
          vec3_t x_cell = cellCentre(full_grid, id_neigh);
          if ((x_cell - x_face)*n_face < 0) { // !!!
            keep_cell[id_neigh] = 1;
          } else {
            if (keep_cell[id_neigh] == 0) {
              keep_cell[id_neigh] = 2;
              seed_cells.append(id_neigh);
            }
          }
        }
      }
    }
  }
  while (seed_cells.size() > 0) {
    QList<vtkIdType> new_seed_cells;
    foreach (vtkIdType id_cell, seed_cells) {
      for (int i = 0; i < full_part.c2cGSize(id_cell); ++i) {
        vtkIdType id_neigh = full_part.c2cGG(id_cell, i);
        if (keep_cell[id_neigh] == 0) {
          new_seed_cells.append(id_neigh);
          keep_cell[id_neigh] = 2;
        }
      }
    }
    seed_cells = new_seed_cells;
  }
  QList<vtkIdType> keep_cells;
  EG_FORALL_CELLS(id_cell, full_grid) {
    if (keep_cell[id_cell] < 2) {
      keep_cells.append(id_cell);
    }
  }
  cout << "keeping " << keep_cells.size() << " cells" << endl;
  MeshPartition keep_part(full_grid);
  keep_part.setCells(keep_cells);
  keep_part.extractToVtkGrid(m_Grid);
  //makeCopy(full_grid, m_Grid);
}

void TetGenOperation::readSettings()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings");
  if(!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaximalEdgeLength;
    in >> m_MinimalEdgeLength;
    in >> m_GrowthFactor;
    in >> m_NodesPerQuarterCircle;
    int num_bcs;
    in >> num_bcs;
    int check_state;
    for (int i = 0; i < num_bcs; ++i) {
      in >> check_state;
    }
    in >> m_2dFeatureResolution;
    in >> m_3dFeatureResolution;
  }

  m_ELSManager.read();
  m_ELSManager.readRules(m_Grid);
}

QString TetGenOperation::qualityText()
{
  using namespace GeometryTools;
  vec3_t x1(1.0/sqrt(3.0), 0, 0);
  vec3_t x2 = rotate(x1, vec3_t(0, 0, 1), deg2rad(60.0));
  vec3_t x3 = rotate(x2, vec3_t(0, 0, 1), deg2rad(60.0));
  vec3_t x4(0, 0, 1.0);
  vec3_t n1 = triNormal(x1, x2, x3);
  vec3_t n2 = triNormal(x1, x4, x2);
  vec3_t n3 = triNormal(x2, x4, x3);
  double alpha = angle(n1, n2);
  alpha = min(alpha, angle(n1, n3));
  alpha = min(alpha, angle(n2, n3));
  QString alpha_txt;
  alpha_txt.setNum(rad2deg(alpha));
  return QString("1.4") + "/" + alpha_txt;
}

void TetGenOperation::tetgen(QString flags, vtkUnstructuredGrid *background_grid)
{
  try {
    tetgenio in, out, background;
    readSettings();
    copyToTetGen(in);
    copyToTetGen(background, background_grid);
    char *char_flags = new char [flags.size() + 1];
    strcpy(char_flags, qPrintable(flags));
    tetrahedralize(char_flags, &in, &out, NULL, &background);
    copyFromTetGen(out);
    delete [] char_flags;
  } catch (int x) {
    switch (x) {
    case 1: // Out of memory.
      EG_ERR_RETURN("TetGen error:  Out of memory");
      break;
    case 2: // Encounter an internal error.
      EG_ERR_RETURN("TetGen error:  internal error");
      break;
    case 3:
      EG_ERR_RETURN("TetGen error:  A self-intersection was detected.");
      break;
    case 4:
      EG_ERR_RETURN("TetGen error:  A very small input feature size was detected.");
      break;
    case 5:
      EG_ERR_RETURN("TetGen error:  Two very close input facets were detected.");
      break;
    case 10:
      EG_ERR_RETURN("TetGen error:  An input error was detected.");
      EG_BUG;
      break;
    default:
      EG_ERR_RETURN("TetGen error:  unknown error");
      break;
    }
  }
}

