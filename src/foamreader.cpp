//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "foamreader.h"
#include "guimainwindow.h"

FoamReader::FoamReader()
{
  EG_TYPENAME;
}

void FoamReader::operate()
{
  try {
    readInputDirectory("Select OpenFOAM case directory");
    setCaseDir(getFileName());
    if (isValid()) {
      buildMaps();
      EG_VTKSP(vtkUnstructuredGrid, ug);
      int num_faces;

      {
        readFile("constant/polyMesh/faces");
        QTextStream f(getBuffer());
        f >> num_faces;
      }
      allocateGrid(ug, num_faces - getFirstBoundaryFace(), numSurfNodes());
      EG_VTKDCC(vtkIntArray, bc, ug, "cell_code");
      EG_VTKDCC(vtkIntArray, orgdir, ug, "cell_orgdir");
      EG_VTKDCC(vtkIntArray, voldir, ug, "cell_voldir");
      EG_VTKDCC(vtkIntArray, curdir, ug, "cell_curdir");
      {
        readFile("constant/polyMesh/points");
        QTextStream f(getBuffer());
        int num_nodes;
        f >> num_nodes;
        for (int i = 0; i < num_nodes; ++i) {
          vec3_t x;
          f >> x[0] >> x[1] >> x[2];
          if (volToSurf(i) != -1) {
            ug->GetPoints()->SetPoint(volToSurf(i), x.data());
          }
        }
      }
      {
        readFile("constant/polyMesh/faces");
        QTextStream f(getBuffer());
        f >> num_faces;
        for (int i = 0; i < num_faces; ++i) {
          vtkIdType N_pts;
          f >> N_pts;
          QVector<vtkIdType> pts(N_pts);
          for (int j = 0; j < N_pts; ++j) {
            f >> pts[j];
          }
          if (i >= getFirstBoundaryFace()) {
            for (int j = 0; j < N_pts; ++j) {
              pts[j] = volToSurf(pts[j]);
            }
            vtkIdType id_cell = 0;
            if (N_pts == 3) {
              id_cell = ug->InsertNextCell(VTK_TRIANGLE, N_pts, pts.data());
            } else if (N_pts == 4) {
              id_cell = ug->InsertNextCell(VTK_QUAD, N_pts, pts.data());
            } else {
              id_cell = ug->InsertNextCell(VTK_POLYGON, N_pts, pts.data());
            }
            bc->SetValue(id_cell, 999);
            orgdir->SetValue(id_cell, 0);
            curdir->SetValue(id_cell, 0);
            voldir->SetValue(id_cell, 0);
          }
        }
      }
      GuiMainWindow::pointer()->clearBCs();
      {
        readFile("constant/polyMesh/boundary");
        QTextStream f(getBuffer());
        int num_patches;
        f >> num_patches;
        for (int i = 1; i <= num_patches; ++i) {
          QString name, dummy, type;
          f >> name >> dummy >> type;
          BoundaryCondition BC(name, type);
          GuiMainWindow::pointer()->addBC(i, BC);
          vtkIdType num_faces, start_face;
          f >> dummy >> num_faces;
          f >> dummy >> start_face;
          for (vtkIdType id_cell = start_face - getFirstBoundaryFace(); id_cell < start_face - getFirstBoundaryFace() + num_faces; ++id_cell) {
            bc->SetValue(id_cell, i);
          }
        }
      }
      makeCopy(ug, m_Grid);
      createBasicFields(m_Grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
      UpdateCellIndex(m_Grid);
    }
  } catch (Error err) {
    err.display();
  }
}
