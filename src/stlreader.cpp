//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "stlreader.h"
#include "correctsurfaceorientation.h"

#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include <vtkSTLReader.h>
#include <vtkCleanPolyData.h>

#include <QFileInfo>
#include "guimainwindow.h"

StlReader::StlReader()
{
  setFormat("STL files(*.stl *.STL)");
};

void StlReader::operate()
{
  QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
  readInputFileName(file_info.completeBaseName() + ".stl");
  if (isValid()) {
    vtkSTLReader *stl = vtkSTLReader::New();
    stl->MergingOn();
    stl->SetFileName(getCFileName());
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
        };
        L = min(L, (x1-x2).abs());
      };
    };
    EG_VTKSP(vtkCleanPolyData, poly_clean);
    poly_clean->ToleranceIsAbsoluteOn();
    poly_clean->SetAbsoluteTolerance(0.5*L);
    poly_clean->SetInput(poly);
    EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, poly2ugrid);
    poly2ugrid->SetInput(poly_clean->GetOutput());
    poly2ugrid->Update();
    
    allocateGrid(grid, poly2ugrid->GetOutput()->GetNumberOfCells(), poly2ugrid->GetOutput()->GetNumberOfPoints());
    for (vtkIdType id_node = 0; id_node < poly2ugrid->GetOutput()->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      poly2ugrid->GetOutput()->GetPoints()->GetPoint(id_node, x.data());
      grid->GetPoints()->SetPoint(id_node, x.data());
    };
    for (vtkIdType id_cell = 0; id_cell < poly2ugrid->GetOutput()->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = poly2ugrid->GetOutput()->GetCellType(id_cell);
      poly2ugrid->GetOutput()->GetCellPoints(id_cell, N_pts, pts);
      grid->InsertNextCell(type_cell, N_pts, pts);
    };
    
    EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
    EG_VTKDCC(vtkIntArray, orgdir, grid, "cell_orgdir");
    EG_VTKDCC(vtkIntArray, voldir, grid, "cell_voldir");
    EG_VTKDCC(vtkIntArray, curdir, grid, "cell_curdir");
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      bc->SetValue(id_cell, 999);
      orgdir->SetValue(id_cell, 0);
      voldir->SetValue(id_cell, 0);
      curdir->SetValue(id_cell, 0);
    };
    CorrectSurfaceOrientation corr_surf;
    corr_surf.setGrid(grid);
    corr_surf();
    
  };
};
