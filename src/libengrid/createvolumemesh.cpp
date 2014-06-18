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
  m_CreateVolumeMesh = false;
  m_BackgroundGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  m_FirstCall = true;
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
      flags += "mR";
    }
    if (preserve_surface) {
      flags += "Y";
    }
    cout << "TetGen pass " << pass << "  flags=" << qPrintable(flags) << endl;
    if (m_FirstCall) {
      tetgen(flags);
      m_FirstCall = false;
      makeCopy(m_Grid, m_BackgroundGrid);
      cout << "background mesh:\n";
      cout << "  nodes: " << m_Grid->GetNumberOfPoints() << endl;
      int N1 = m_Grid->GetNumberOfCells();
      int N2 = numVolumeCells();
      cout << "  faces: " << N1 - N2 << endl;
      cout << "  cells: " << N2 << endl;
    } else {
      tetgen(flags, m_BackgroundGrid);
    }
    N2 = numVolumeCells();
    cout << N2 << endl;
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

void CreateVolumeMesh::setVolumeMeshOn()
{
  m_CreateVolumeMesh = true;
}

void CreateVolumeMesh::reduceSurface(QSet<int> boundary_codes)
{
  RemovePoints remove_points;
  MeshPartition part;
  part.setGrid(m_Grid);
  part.setAllCells();
  remove_points.setMeshPartition(part);
  remove_points.setBoundaryCodes(boundary_codes);
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


void CreateVolumeMesh::operate()
{
  readSettings();
  m_Part.trackGrid(m_Grid);

  if (m_CreateBoundaryLayer) {
    cout << "A" << endl;
    createTetMesh(2, false);
    CreateBoundaryLayerShell blayer;
    blayer.setGrid(m_Grid);
    blayer.setAllCells();
    blayer();
    if (m_CreateVolumeMesh) {
      cout << "B" << endl;
      createTetMesh(1, true);
      vtkUnstructuredGrid *prismatic_grid = blayer.getPrismaticGrid();
      MeshPartition prismatic_part(prismatic_grid, true);
      QVector<vtkIdType> shell_cells;
      getSurfaceCells(blayer.getBoundaryLayerCodes(), shell_cells, m_Grid);
      DeleteCells delete_cells;
      delete_cells.setGrid(m_Grid);
      delete_cells.setCellsToDelete(shell_cells);
      delete_cells();
      m_Part.addPartition(prismatic_part);
    }
  } else if (m_CreateVolumeMesh) {
    cout << "C" << endl;
    createTetMesh(3, false);
  }
}

