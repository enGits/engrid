//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#include "su2writer.h"
#include "guimainwindow.h"

Su2Writer::Su2Writer()
{
  setFormat("SU2 mesh files(*.su2)");
}

void Su2Writer::writeHeader()
{
  QTextStream f(m_File);
  f << "%\n";
  f << "% Problem dimension\n";
  f << "%\n";
  f << "NDIME= 3\n";
}

void Su2Writer::writeElements()
{
  QTextStream f(m_File);
  f << "%\n";
  f << "% Inner element connectivity\n";
  f << "%\n";
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      ++N;
    }
  }
  f << "NELEM= " << N << "\n";
  int i = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      f << m_Grid->GetCellType(id_cell);
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int j = 0; j < N_pts; ++j) {
        f << " " << pts[j];
      }
      f << " " << i << "\n";
      ++i;
    }
  }
}

void Su2Writer::writeNodes()
{
  QTextStream f(m_File);
  f << "%\n";
  f << "% Node coordinates\n";
  f << "%\n";
  f << "NPOIN= " << m_Grid->GetNumberOfPoints() << "\n";
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << " " << id_node << "\n";
  }
}

void Su2Writer::writeBoundaries()
{
  QTextStream f(m_File);
  f << "%\n";
  f << "% Boundary elements\n";
  f << "%\n";
  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  f << "NMARK= " << bcs.size() << "\n";
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (int bc, bcs) {
    BoundaryCondition BC = GuiMainWindow::pointer()->getBC(bc);
    f << "MARKER_TAG= " << BC.getName() << "\n";
    QList<vtkIdType> faces;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (cell_code->GetValue(id_cell) == bc) {
          faces.append(id_cell);
        }
      }
    }
    f << "MARKER_ELEMS= " << faces.size() << "\n";
    foreach (vtkIdType id_cell, faces) {
      f << m_Grid->GetCellType(id_cell);
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int j = 0; j < N_pts; ++j) {
        f << " " << pts[j];
      }
      f << "\n";
    }
  }
}

void Su2Writer::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".su2");
    if (isValid()) {
      m_File = new QFile(getFileName());
      m_File->open(QIODevice::WriteOnly | QIODevice::Text);
      writeHeader();
      writeElements();
      writeNodes();
      writeBoundaries();
      delete m_File;
    }
  } catch (Error err) {
    err.display();
  }
}

