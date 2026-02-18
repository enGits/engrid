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

#include "createvolumemesh.h"
#include "deletetetras.h"
#include "deletecells.h"
#include "engrid.h"
#include "guimainwindow.h"
#include "meshpartition.h"
#include "updatedesiredmeshdensity.h"
#include "createboundarylayershell.h"

#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

CreateVolumeMesh::CreateVolumeMesh()
{
  EG_TYPENAME;
  m_CreateBoundaryLayer = false;
  m_BackgroundGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_FirstCall = true;
  m_Debug = false;
}

int CreateVolumeMesh::numVolumeCells()
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      ++N;
    }
  }
  return N;
}

void CreateVolumeMesh::createTetMesh(int max_num_passes, bool preserve_surface)
{
  double a = m_MaximalEdgeLength;
  double V = a*a*a/(6*sqrt(2.0));
  int N1 = 0;
  int N2 = numVolumeCells();
  int pass = 1;
  QString V_txt;
  V_txt.setNum(V);
  bool done = false;
  while (!done) {
    N1 = N2;
    QString flags;
    QString q_txt = "1.2/0";
    if (pass > 1) {
      q_txt = "1.2/0";
    }
    flags = QString("pq") + q_txt + "a" + V_txt;
    if (!m_FirstCall) {
      flags += "m";
    }
    if (preserve_surface) {
      flags += "Y";
    }
    cout << "TetGen pass " << pass << "  flags=" << qPrintable(flags) << endl;
    if (m_FirstCall) {
      m_FirstCall = false;
      tetgen(flags);
    } else {
      makeCopy(m_Grid, m_BackgroundGrid);
      tetgen(flags, m_BackgroundGrid);
    }

    N2 = numVolumeCells();
    ++pass;
    if (fabs(double(N2-N1)/N1) < 0.05 || pass > max_num_passes) {
      done = true;
    }
  }
}

void CreateVolumeMesh::setBoundaryLayerOn()
{
  m_CreateBoundaryLayer = true;
}

void CreateVolumeMesh::triangulateBoundaryLayerInterface(vtkUnstructuredGrid *shell_grid)
{
  EG_VTKDCC(vtkIntArray, cell_code, shell_grid, "cell_code");

  vtkIdType num_tri_cells = 0;
  for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
    if (!isSurface(id_cell, shell_grid)) {
      continue;
    }
    EG_GET_CELL(id_cell, shell_grid);
    if (num_pts >= 3) {
      if (num_pts == 3) {
        ++num_tri_cells;
      } else {
        num_tri_cells += num_pts - 2;
      }
    }
  }

  EG_VTKSP(vtkUnstructuredGrid, tri_grid);
  allocateGrid(tri_grid, num_tri_cells, shell_grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < shell_grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    shell_grid->GetPoint(id_node, x.data());
    tri_grid->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(shell_grid, id_node, tri_grid, id_node);
  }

  EG_VTKDCC(vtkIntArray, tri_cell_code, tri_grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
    if (!isSurface(id_cell, shell_grid)) {
      continue;
    }
    EG_GET_CELL(id_cell, shell_grid);
    if (num_pts < 3) {
      continue;
    }
    if (num_pts == 3) {
      vtkIdType tri_pts[3] = { pts[0], pts[1], pts[2] };
      vtkIdType id_new = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
      copyCellData(shell_grid, id_cell, tri_grid, id_new);
      tri_cell_code->SetValue(id_new, cell_code->GetValue(id_cell));
    } else {
      QVector<vtkIdType> poly;
      poly.resize(num_pts);
      for (int i = 0; i < num_pts; ++i) {
        poly[i] = pts[i];
      }
      QVector<QVector<vtkIdType> > triangles;
      triangulatePolygon(shell_grid, poly, triangles);
      for (int i = 0; i < triangles.size(); ++i) {
        vtkIdType tri_pts[3] = { triangles[i][0], triangles[i][1], triangles[i][2] };
        vtkIdType id_new = tri_grid->InsertNextCell(VTK_TRIANGLE, 3, tri_pts);
        copyCellData(shell_grid, id_cell, tri_grid, id_new);
        tri_cell_code->SetValue(id_new, cell_code->GetValue(id_cell));
      }
    }
  }

  makeCopy(tri_grid, shell_grid);
}

void CreateVolumeMesh::operate()
{
  readSettings();
  m_Part.trackGrid(m_Grid);

  if (m_CreateBoundaryLayer) {
    CreateBoundaryLayerShell blayer;
    blayer.setGrid(m_Grid);
    blayer.setAllCells();
    blayer();
    if (!m_Debug) {
      if (blayer.success()) {
        //createTetMesh(2, true);
        //
        vtkUnstructuredGrid *prismatic_grid = blayer.getPrismaticGrid();
        MeshPartition prismatic_part(prismatic_grid, true);
        QVector<vtkIdType> shell_cells;
        getSurfaceCells(blayer.getBoundaryLayerCodes(), shell_cells, m_Grid);
        MeshPartition shell_part(m_Grid);
        shell_part.setCells(shell_cells);
        EG_VTKSP(vtkUnstructuredGrid, shell_grid);
        shell_part.extractToVtkGrid(shell_grid);
        if (blayer.usePolyhedralBoundaryLayer()) {
          triangulateBoundaryLayerInterface(shell_grid);
        }
        //
        auto new_bc = GuiMainWindow::pointer()->getBC(BoundaryCondition("transfer", "patch"));
        {
          EG_VTKDCC(vtkIntArray, cell_code, shell_grid, "cell_code");
          for (vtkIdType id_cell = 0; id_cell < shell_grid->GetNumberOfCells(); ++id_cell) {
            //
            // set boundary code for shell cells (outside end of boundary layer mesh)
            //
            cell_code->SetValue(id_cell, new_bc.getCode());
            //
            // invert the direction of the normals of the shell cells
            //
            shell_grid->GetCells()->ReverseCellAtId(id_cell);
          }
        }
        //
        DeleteCells delete_cells;
        MeshPartition all(m_Grid, true);
        delete_cells.setGrid(m_Grid);
        delete_cells.setCellsToDelete(all.getCells());
        delete_cells();
        m_Part.setAllCells();
        m_Part.addPartition(prismatic_part);
        //
        MeshPartition new_shell_part(shell_grid, true);
        m_Part.addPartition(new_shell_part);
        //
      } else {
        cout << "An error occurred while creating the prismatic layers!" << endl;
      }
    }
  } else {
    createTetMesh(2, true);
  }
  m_Grid->Modified();
}
