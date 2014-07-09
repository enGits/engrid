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

#include "checkforoverlap.h"
#include "guimainwindow.h"
#include "timer.h"
#include "facefinder.h"

CheckForOverlap::CheckForOverlap()
{
  EG_TYPENAME;
}

void CheckForOverlap::operate()
{
  cout << "testing mesh for overlapping faces ..." << endl;
  Timer timer;
  m_BoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();
  QVector<bool> is_overlap(m_Grid->GetNumberOfCells(), false);
  QVector<vtkIdType> faces;
  getAllSurfaceCells(faces, m_Grid);
  setAllSurfaceCells();
  int bc_max = 0;  
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");

  FaceFinder find;
  find.setMaxNumFaces(100);
  find.setGrid(m_Grid);

  int N_searches = 0;
  int N_buckets = 0;
  int N_faces = 0;
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {    
    vtkIdType id_face1 = faces[i_faces];
    bc_max = max(bc_max, cell_code->GetValue(id_face1));
    QVector<vec3_t> x1;
    vec3_t xc(0,0,0);
    QSet<vtkIdType> neighbours;
    {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_face1, N_pts, pts);
      x1.resize(N_pts + 1);
      for (int i = 0; i < N_pts; ++i) {
        m_Grid->GetPoint(pts[i], x1[i].data());
        xc += x1[i];
        for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
          neighbours.insert(m_Part.n2cGG(pts[i], j));
        }
      }
      xc *= 1.0/N_pts;
      x1[N_pts] = x1[0];
    }
    QVector<vtkIdType> close_faces;
    find.getCloseFaces(xc, close_faces);
    if (close_faces.size() == 0) {
      EG_BUG;
    }
    ++N_buckets;
    N_faces += close_faces.size();
    foreach (vtkIdType id_face2, close_faces) {
      if (!neighbours.contains(id_face2)) {
        ++N_searches;
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_face2, N_pts, pts);
        QVector<vec3_t> x2(N_pts);
        for (int i = 0; i < N_pts; ++i) {
          m_Grid->GetPoint(pts[i], x2[i].data());
        }
        for (int i = 0; i < x1.size() - 1; ++i) {
          vec3_t x, r;
          if (intersectEdgeAndTriangle(x2[0], x2[1], x2[2], x1[i], x1[i+1], x, r, 1e-6)) {
            is_overlap[id_face1] = true;
            is_overlap[id_face2] = true;
          }
        }
      }
    }
    if (timer()) {
      cout << "  " << i_faces + 1 << " of " << faces.size() << "  " << N_searches << " searches" << endl;
      cout << N_faces/N_buckets << " faces/bucket" << endl;
      N_searches = 0;
    }
  }
  int N = 0;
  foreach (vtkIdType id_face, faces) {
    if (is_overlap[id_face]) {
      cell_code->SetValue(id_face, bc_max + 1);
      ++N;
    }
  }
  if (N > 0) {
    GuiMainWindow::pointer()->addBC(bc_max + 1, BoundaryCondition("overlapping_faces", "patch"));
    GuiMainWindow::pointer()->updateBoundaryCodes(true);
  }
  cout << N << " overlapping or close faces found" << endl;
  cout << N_faces/N_buckets << " faces/bucket" << endl;
}
