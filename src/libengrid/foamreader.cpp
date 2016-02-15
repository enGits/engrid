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
      buildFoamMaps();
      EG_VTKSP(vtkUnstructuredGrid, ug);
      int num_faces;

      {
        readFile("constant/polyMesh/faces");
        QTextStream f(getBuffer());
        f >> num_faces;
      }

      // count faces
      int num_triangles = 0;
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
            if (N_pts < 3) {
              EG_BUG;
            }
            num_triangles += (N_pts - 2);
          }
        }
      }

      allocateGrid(ug, num_triangles, numSurfNodes());
      EG_VTKDCC(vtkIntArray, bc, ug, "cell_code");
      EG_VTKDCC(vtkIntArray, orgdir, ug, "cell_orgdir");
      EG_VTKDCC(vtkIntArray, voldir, ug, "cell_voldir");
      EG_VTKDCC(vtkIntArray, curdir, ug, "cell_curdir");
      QVector<vtkIdType> face2id(num_faces, -1);
      QVector<int> num_tris(num_faces, 0);

      // insert points
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

      qWarning() << ug->GetNumberOfPoints() << " nodes";
      qWarning() << num_faces - getFirstBoundaryFace() << " faces";
      qWarning() << num_triangles << " triangles";
      qWarning() << num_tris.size();

      // insert faces
      {
        readFile("constant/polyMesh/faces");
        QTextStream f(getBuffer());
        f >> num_faces;
        for (int i = 0; i < num_faces; ++i) {
          vtkIdType N_pts;
          f >> N_pts;
          QVector<vtkIdType> pts(N_pts+2);
          for (int j = 0; j < N_pts; ++j) {
            f >> pts[j];
          }
          if (i >= getFirstBoundaryFace()) {
            for (int j = 0; j < N_pts; ++j) {
              pts[j] = volToSurf(pts[j]);
            }
            pts[N_pts]   = pts[0];
            pts[N_pts+1] = pts[1];
            QVector<vec3_t> x(N_pts+2);
            vec3_t xc(0,0,0);
            for (int j = 0; j < N_pts; ++j) {
              ug->GetPoint(pts[j], x[j].data());
              xc += x[j];
            }
            x[N_pts]   = x[0];
            x[N_pts+1] = x[1];
            xc *= 1.0/N_pts;
            vec3_t n(0,0,0);
            for (int j = 0; j < N_pts; ++j) {
              vec3_t u = x[j] - xc;
              vec3_t v = x[j+1] - xc;
              n += u.cross(v);
            }
            n.normalise();
            int j0 = 0;
            double scal0 = 1e99;
            for (int j = 1; j <= N_pts; ++j) {
              vec3_t u = x[j]   - x[j-1];
              vec3_t v = x[j+1] - x[j];
              u.normalise();
              v.normalise();
              vec3_t nj = u.cross(v);
              double scal = nj*n;
              if (scal < scal0) {
                j0 = j % N_pts;
                scal0 = scal;
              }
            }
            QList<int> jn;
            for (int j = 0; j < j0 - 1; ++j) {
              jn.append(j);
            }
            if (j0 > 0) {
              for (int j = j0 + 1; j < N_pts; ++j) {
                jn.append(j);
              }
            } else {
              for (int j = j0 + 1; j < N_pts - 1; ++j) {
                jn.append(j);
              }
            }
            foreach (int j, jn) {
              vtkIdType p[3];
              p[0] = pts[j];
              p[1] = pts[j+1];
              p[2] = pts[j0];
              vtkIdType id_cell = ug->InsertNextCell(VTK_TRIANGLE, 3, p);
              bc->SetValue(id_cell, 999);
              orgdir->SetValue(id_cell, 0);
              curdir->SetValue(id_cell, 0);
              voldir->SetValue(id_cell, 0);
              if (num_tris[i] == 0) {
                face2id[i] = id_cell;
              }
              ++num_tris[i];
            }
          }
        }
      }
      GuiMainWindow::pointer()->clearBCs();
      {
        int N1 = 0;
        int N2 = 0;
        readFile("constant/polyMesh/boundary");
        QTextStream f(getBuffer());
        int num_patches;
        f >> num_patches;
        for (int i = 1; i <= num_patches; ++i) {
          QString name, dummy, type;
          f >> name >> dummy >> type;
          BoundaryCondition BC(name, type, i);
          GuiMainWindow::pointer()->setBC(i, BC);
          vtkIdType num_faces, start_face;
          f >> dummy >> num_faces;
          f >> dummy >> start_face;
          for (vtkIdType i_face = start_face; i_face < start_face + num_faces; ++i_face) {
            if (face2id[i_face] < 0) {
              EG_BUG;
            }
            ++N1;
            N2 += num_tris[i_face];
            for (int j = 0; j < num_tris[i_face]; ++j) {
              vtkIdType id_cell = face2id[i_face] + j;
              bc->SetValue(id_cell, i);
            }
          }
          qWarning() << i << ',' << N1 << ',' << N2;
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
