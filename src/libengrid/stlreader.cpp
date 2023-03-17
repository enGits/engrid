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
#include "stlreader.h"
#include "correctsurfaceorientation.h"

#include "engrid.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include <vtkIdList.h>
#include <vtkSTLReader.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>

#include <QFileInfo>
#include <QInputDialog>
#include <vtkSmartPointer.h>

#include "guimainwindow.h"
#include "fixcadgeometry.h"

StlReader::StlReader()
{
  setFormat("STL files(*.stl *.STL)");
  m_Tolerance = -1;
  m_FileNameSet = false;
  m_MaxNumCleanIter = 20;
}

void StlReader::setFileName(QString file_name)
{
  m_FileName = file_name;
  cout << qPrintable(m_FileName) << endl;
  m_FileNameSet = true;
}

void StlReader::operate()
{
  QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
  QString file_name;
  if (m_FileNameSet) {
    file_name = m_FileName;
  } else {
    readInputFileName(file_info.completeBaseName() + ".stl");
    if (isValid()) {
      file_name = getFileName();
    } else {
      return;
    }
  }

  EG_VTKSP(vtkSTLReader, stl);
  stl->MergingOn();
  stl->SetFileName(qPrintable(file_name));
  stl->Update();
  EG_VTKSP(vtkPolyData, poly);
  poly->DeepCopy(stl->GetOutput());
  poly->BuildCells();
  double L = 1e99;
  for (vtkIdType cellId = 0; cellId < poly->GetNumberOfCells(); ++cellId) {
    EG_GET_CELL(cellId, poly);
    for (int i = 0; i < num_pts; ++i) {
      vec3_t x1, x2;
      poly->GetPoints()->GetPoint(pts[i], x1.data());
      if (i == num_pts - 1) {
        poly->GetPoints()->GetPoint(pts[0], x2.data());
      } else {
        poly->GetPoints()->GetPoint(pts[i+1], x2.data());
      }
      L = min(L, (x1-x2).abs());
    }
  }
  if (m_Tolerance < 0) {
    m_Tolerance = QInputDialog::getText(NULL, "enter STL tolerance", "tolerance", QLineEdit::Normal, "1e-10").toDouble();
  }
  cout << "cleaning STL geometry:" << endl;
  EG_VTKSP(vtkCleanPolyData, poly_clean);
  EG_VTKSP(vtkFeatureEdges, topo_check);
  double bounds[6];
  poly->GetBounds(bounds);
  poly_clean->ToleranceIsAbsoluteOn();
  poly_clean->ConvertLinesToPointsOn();
  poly_clean->ConvertPolysToLinesOn();
  poly_clean->SetInputData(poly);
  topo_check->SetInputConnection(poly_clean->GetOutputPort());
  topo_check->BoundaryEdgesOn();
  topo_check->ManifoldEdgesOff();
  topo_check->FeatureEdgesOff();
  topo_check->NonManifoldEdgesOn();
  bool check_passed;
  int count = 0;
  do {
    ++count;
    cout << "  tolerance = " << m_Tolerance << endl;
    poly_clean->SetAbsoluteTolerance(m_Tolerance);
    topo_check->Update();
    m_Tolerance *= 1.5;
    check_passed = topo_check->GetOutput()->GetNumberOfPoints() == 0;
  } while (m_Tolerance < 1 && !check_passed && count < m_MaxNumCleanIter);
  if (check_passed) {
    cout << "The STL geometry seems to be clean." << endl;
  } else {
    cout << "The STL geometry could not be cleaned." << endl;
  }
  // with a tolerance of " << 0.5*L << endl;
  EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, poly2ugrid);
  poly2ugrid->SetInputConnection(poly_clean->GetOutputPort());
  poly2ugrid->Update();

  allocateGrid(m_Grid, poly2ugrid->GetOutput()->GetNumberOfCells(), poly2ugrid->GetOutput()->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < poly2ugrid->GetOutput()->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    poly2ugrid->GetOutput()->GetPoints()->GetPoint(id_node, x.data());
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
  for (vtkIdType id_cell = 0; id_cell < poly2ugrid->GetOutput()->GetNumberOfCells(); ++id_cell) {
    vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
    vtkIdType type_cell = poly2ugrid->GetOutput()->GetCellType(id_cell);
    poly2ugrid->GetOutput()->GetCellPoints(id_cell, pts);
    m_Grid->InsertNextCell(type_cell, pts);
  }

  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  EG_VTKDCC(vtkIntArray, orgdir, m_Grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, voldir, m_Grid, "cell_voldir");
  EG_VTKDCC(vtkIntArray, curdir, m_Grid, "cell_curdir");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    bc->SetValue(id_cell, 1);
    orgdir->SetValue(id_cell, 0);
    voldir->SetValue(id_cell, 0);
    curdir->SetValue(id_cell, 0);
  }
  if (check_passed) {
    CorrectSurfaceOrientation corr_surf;
    corr_surf.setGrid(m_Grid);
    corr_surf();
    FixCadGeometry cad_fix;
    cad_fix.setGrid(m_Grid);
    cad_fix();
  }
}
