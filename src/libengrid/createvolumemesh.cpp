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
#include "guimainwindow.h"
#include "updatedesiredmeshdensity.h"
#include "createboundarylayershell.h"

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
    QString q_txt = "1.2/18";
    if (pass > 1) {
      q_txt = "1.2/18";
    }
    flags = QString("pq") + q_txt + "a" + V_txt + "S1000000";
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
        createTetMesh(2, true);
        vtkUnstructuredGrid *prismatic_grid = blayer.getPrismaticGrid();
        MeshPartition prismatic_part(prismatic_grid, true);
        QVector<vtkIdType> shell_cells;
        getSurfaceCells(blayer.getBoundaryLayerCodes(), shell_cells, m_Grid);
        DeleteCells delete_cells;
        delete_cells.setGrid(m_Grid);
        delete_cells.setCellsToDelete(shell_cells);
        delete_cells();
        m_Part.addPartition(prismatic_part);
      } else {
        cout << "An error ocuured while creating the prismatic layers!" << endl;
      }
    }
  } else {
    createTetMesh(2, true);
  }
}

