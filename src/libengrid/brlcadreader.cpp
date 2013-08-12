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

#ifdef BRLCAD_SUPPORT

#include "brlcadreader.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include "setboundarycode.h"
#include "facefinder.h"
#include "guimainwindow.h"
#include "guibrlcadimportdialogue.h"
#include "createcadtesselation.h"
#include "stlreader.h"
#include "brlcadinterface.h"

#include <vtkSTLReader.h>

BrlcadReader::BrlcadReader()
{
  EG_TYPENAME;
  setFormat("BRL-CAD database files(*.g)");
  setExtension(".g");
}

void BrlcadReader::processStlFile(QString file_name, bool append_to_list)
{
  vtkSTLReader *stl = vtkSTLReader::New();
  stl->MergingOn();
  stl->SetFileName(file_name.toAscii().data());
  stl->Update();
  EG_VTKSP(vtkPolyData, poly);
  poly->DeepCopy(stl->GetOutput());
  poly->BuildCells();
  double L = 1e99;
  for (vtkIdType cellId = 0; cellId < poly->GetNumberOfCells(); ++cellId) {
    vtkIdType *pts, Npts;
    poly->GetCellPoints(cellId, Npts, pts);
    for (int i = 0; i < Npts; ++i) {
      vec3_t x1, x2;
      poly->GetPoints()->GetPoint(pts[i], x1.data());
      if (i == Npts - 1) {
        poly->GetPoints()->GetPoint(pts[0], x2.data());
      } else {
        poly->GetPoints()->GetPoint(pts[i+1], x2.data());
      }
      L = min(L, (x1-x2).abs());
    }
  }
  EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, poly2ugrid);
  poly2ugrid->SetInput(poly);
  poly2ugrid->Update();

  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, poly2ugrid->GetOutput()->GetNumberOfCells(), poly2ugrid->GetOutput()->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < poly2ugrid->GetOutput()->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    poly2ugrid->GetOutput()->GetPoints()->GetPoint(id_node, x.data());
    grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (vtkIdType id_cell = 0; id_cell < poly2ugrid->GetOutput()->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = poly2ugrid->GetOutput()->GetCellType(id_cell);
    poly2ugrid->GetOutput()->GetCellPoints(id_cell, N_pts, pts);
    grid->InsertNextCell(type_cell, N_pts, pts);
  }

  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCC(vtkIntArray, orgdir,    grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, voldir,    grid, "cell_voldir");
  EG_VTKDCC(vtkIntArray, curdir,    grid, "cell_curdir");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cell_code->SetValue(id_cell, 9999);
    orgdir->SetValue(id_cell, 0);
    voldir->SetValue(id_cell, 0);
    curdir->SetValue(id_cell, 0);
  }

  if (append_to_list) {
    SetBoundaryCode set_bc;
    set_bc.setGrid(grid);
    set_bc.setAllSurfaceCells();
    int bc_max = 1;
    bool done = false;
    do {
      vtkIdType id_start = -1;
      for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
        if (cell_code->GetValue(id_cell) == 9999) {
          id_start = id_cell;
          break;
        }
      }
      if (id_start == -1) {
        done = true;
      } else {
        set_bc.setFeatureAngle(20.0);
        set_bc.setBC(bc_max);
        set_bc.setProcessAll(true);
        set_bc.setSelectAllVisible(false);
        set_bc.setOnlyPickedCell(false);
        set_bc.setOnlyPickedCellAndNeighbours(false);
        set_bc.setStart(id_start);
        set_bc();
        ++bc_max;
      }
    } while (!done);
    cout << "file: " << qPrintable(file_name) << endl;
    for (int bc = 1; bc < bc_max; ++bc) {
      QList<vtkIdType> cells;
      for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
        if (cell_code->GetValue(id_cell) == bc) {
          cells.append(id_cell);
        }
      }
      vtkUnstructuredGrid *new_grid = vtkUnstructuredGrid::New();
      QString bc_txt;
      bc_txt.setNum(bc);
      if (bc_max >= 10) {
        bc_txt = bc_txt.rightJustified(2, '0');
      }
      QFileInfo file_info(file_name);
      m_BCNames[new_grid] = file_info.baseName() + "." + bc_txt;
      cout << "  " << qPrintable(file_info.baseName() + "." + bc_txt) << endl;
      makeCopy(grid, new_grid, cells);
      m_Grids.append(new_grid);
    }
  } else {
    makeCopy(grid, m_Grid);
  }
}

void BrlcadReader::findBoundaryCodes()
{
  int num_grids = m_Grids.size();
  QVector<vtkUnstructuredGrid*> grids(num_grids);
  qCopy(m_Grids.begin(), m_Grids.end(), grids.begin());
  QVector<FaceFinder> finders(num_grids);
  for (int i_grid = 0; i_grid < num_grids; ++i_grid) {
    finders[i_grid].setGrid(grids[i_grid]);
  }
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  bool has_errors = false;
  QVector<bool> bc_exists(num_grids, false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    int best_grid = -1;
    double L_min = 1e99;
    vec3_t x1 = cellCentre(m_Grid, id_cell);
    for (int i_grid = 0; i_grid < num_grids; ++i_grid) {
      double L = 1e99;
      vtkIdType id_closest = finders[i_grid].getClosestFace(x1, L);
      vec3_t x2(-999,-999,-999);
      if (id_closest != -1) x2 = cellCentre(m_Grids[i_grid], id_closest);
      if (id_closest != -1) {
        if (L < L_min) {
          best_grid = i_grid;
          L_min = L;
        }
      }
    }
    if (best_grid == -1) {
      has_errors = true;
      cell_code->SetValue(id_cell, 9999);
    } else {
      bc_exists[best_grid] = true;
      cell_code->SetValue(id_cell, best_grid + 1);
    }
  }
  GuiMainWindow::pointer()->clearBCs();
  int bc_max = 1;
  QVector<int> bc_map(num_grids+1,9999);
  m_BC2GridIndex.clear();
  for (int i_grid = 0; i_grid < num_grids; ++i_grid) {
    if (bc_exists[i_grid]) {
      bc_map[i_grid+1] = bc_max;
      GuiMainWindow::pointer()->addBC(bc_max, BoundaryCondition(m_BCNames[grids[i_grid]], "patch"));
      m_BC2GridIndex[bc_max] = i_grid;
      ++bc_max;
    }
  }
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (cell_code->GetValue(id_cell) != 9999) {
      cell_code->SetValue(id_cell, bc_map[cell_code->GetValue(id_cell)]);
    }
  }
  if (has_errors) {
    GuiMainWindow::pointer()->addBC(9999, BoundaryCondition("error-faces", "patch"));
  }
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
}

void BrlcadReader::createBackgroundGeometry()
{
  // make a backup of the existing grid
  EG_VTKSP(vtkUnstructuredGrid, backup_grid);
  makeCopy(m_Grid, backup_grid);

  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();

  // count total number of nodes and faces for background geometry
  int num_nodes = 0;
  int num_faces = 0;
  foreach (int bc, bcs) {
    int i_grid = m_BC2GridIndex[bc];
    num_nodes += m_Grids[i_grid]->GetNumberOfPoints();
    num_faces += m_Grids[i_grid]->GetNumberOfCells();
  }
  allocateGrid(m_Grid, num_faces, num_nodes);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  // copy STL grids into background grid
  vtkIdType id_node_max = 0;
  foreach (int bc, bcs) {
    int i_grid = m_BC2GridIndex[bc];

    // copy nodes
    QVector<int> local2global(m_Grids[i_grid]->GetNumberOfPoints());
    for (vtkIdType id_node = 0; id_node < m_Grids[i_grid]->GetNumberOfPoints(); ++id_node) {
      local2global[id_node] = id_node_max;
      vec3_t x;
      m_Grids[i_grid]->GetPoint(id_node, x.data());
      m_Grid->GetPoints()->SetPoint(id_node_max, x.data());
      ++id_node_max;
    }

    // copy faces
    for (vtkIdType id_face = 0; id_face < m_Grids[i_grid]->GetNumberOfCells(); ++id_face) {
      vtkIdType *pts, N_pts;
      m_Grids[i_grid]->GetCellPoints(id_face, N_pts, pts);
      QVector<vtkIdType> new_pts(N_pts);
      for (int i = 0; i < N_pts; ++i) {
        new_pts[i] = local2global[pts[i]];
      }
      vtkIdType id_new_face = m_Grid->InsertNextCell(m_Grids[i_grid]->GetCellType(id_face), N_pts, new_pts.data());
      cell_code->SetValue(id_new_face, bc);
    }
  }

  // store background grid for surface projection
  GuiMainWindow::pointer()->storeCadInterfaces();

  // restore the initial grid from the backup copy
  makeCopy(backup_grid, m_Grid);
}

void BrlcadReader::operateOld()
{
  readInputDirectory("Select BRL-CAD export directory");
  if (isValid()) {
    QDir dir(getFileName());
    QStringList filters;
    filters << "*.s.stl";
    QStringList file_names = dir.entryList(filters);
    foreach (QString file_name, file_names) {
      processStlFile(dir.path() + "/"+ file_name);
    }
    processStlFile(dir.path() + "/volume.stl", false);
    findBoundaryCodes();
    createBackgroundGeometry();
  }
}

void BrlcadReader::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readInputFileName(file_info.completeBaseName() + ".g");
    if (isValid()) {
      GuiBrlCadImportDialogue dlg;
      dlg.prepare(getFileName());
      if (dlg.exec()) {
        BrlCadInterface* brlcad_interface = new BrlCadInterface(getFileName(), dlg.selectedObject());
        if (dlg.useStlFile()) {
          StlReader stl;
          stl.setFileName(dlg.stlFileName());
          stl();
        } else if (dlg.hasSelectedObject()) {
          CreateCadTesselation tess(brlcad_interface);
          tess.setScanMemory(dlg.scanMemory());
          tess.setPreservationType(dlg.preservationType());
          tess.setSmoothingIterations(dlg.smoothingIterations());
          tess.setSmallestFeatureSize(dlg.smallestFeatureSize());
          tess.setSmallestResolution(dlg.smallestResolution());
          tess.setTargetReduction(dlg.reduction());
          tess();
        }
        GuiMainWindow::pointer()->setUniversalCadInterface(brlcad_interface);
      }
    }
  } catch (Error err) {
    err.display();
  }
}

#endif // BRLCAD_SUPPORT
