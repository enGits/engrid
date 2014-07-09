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

#include "dolfynwriter.h"
#include "guimainwindow.h"

DolfynWriter::DolfynWriter()
{
  setFormat("Dolfyn mesh files(*.vrt)");
}

void DolfynWriter::writeVertices()
{
  QTextStream f(m_VrtFile);
  QString str;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    str.sprintf("%9d      %16.9E%16.9E%16.9E\n", id_node + 1, x[0], x[1], x[2]);
    f << str;
  }
}

void DolfynWriter::writeElements()
{
  int elid = 1;
  QTextStream f(m_CelFile);
  QString str;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType cellId = 0; cellId < m_Grid->GetNumberOfCells(); ++cellId) {
    vtkIdType  Npts;
    vtkIdType *pts;
    m_Grid->GetCellPoints(cellId, Npts, pts);
    if (m_Grid->GetCellType(cellId) == VTK_HEXAHEDRON) {
        str.sprintf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %4d %4d\n",
                    elid, 
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[3] + 1,
                    pts[4] + 1, pts[5] + 1, pts[6] + 1, pts[7] + 1,
                    1, 1);
    }  else if (m_Grid->GetCellType(cellId) == VTK_TETRA) {
        str.sprintf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %4d %4d\n",
                    elid, 
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[2] + 1,
                    pts[3] + 1, pts[3] + 1, pts[3] + 1, pts[3] + 1,
                    1, 1);
    } else if (m_Grid->GetCellType(cellId) == VTK_PYRAMID) {
        str.sprintf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %4d %4d\n",
                    elid, 
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[3] + 1,
                    pts[4] + 1, pts[4] + 1, pts[4] + 1, pts[4] + 1,
                    1, 1);
    } else if (m_Grid->GetCellType(cellId) == VTK_WEDGE) {
        /*
        str.sprintf("%9d      %9d%9d%9d%9d%9d%9d%9d%9d    %5d%5d\n",
                    elid, 
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[2] + 1,
                    pts[3] + 1, pts[4] + 1, pts[5] + 1, pts[5] + 1,
                    1, 1); */
        str.sprintf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %4d %4d\n",
                    elid, 
                    pts[3] + 1, pts[4] + 1, pts[5] + 1, pts[5] + 1,
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[2] + 1,
                    1, 1);
    } else {
        //f << "Skipped!\n"
        continue;
    }
    f << str;
    elid++;
  }
}

void DolfynWriter::writeBoundaries()
{
  int bndid = 1, bcid = 1;
  QTextStream f(m_BndFile);
  QString str;
  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  foreach (int bc, bcs) {
    BoundaryCondition BC = GuiMainWindow::pointer()->getBC(bc);
    QString bc_name = BC.getName();
    QList<vtkIdType> faces;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (isSurface(id_cell, m_Grid)) {
        if (cell_code->GetValue(id_cell) == bc) {
          faces.append(id_cell);
        }
      }
    }
    foreach (vtkIdType id_cell, faces) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      if (m_Grid->GetCellType(id_cell) == VTK_TRIANGLE) {
        str.sprintf("%8d %8d %8d %8d %9d %4d %4d",
                    bndid,
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[2] + 1, 
                    bcid, 0);
      } else if (m_Grid->GetCellType(id_cell) == VTK_QUAD) {
        str.sprintf("%8d %8d %8d %8d %9d %4d %4d",
                    bndid,
                    pts[0] + 1, pts[1] + 1, pts[2] + 1, pts[3] + 1, 
                    bcid, 0);
      } else {
        //f << "Skipped!\n";
        continue;
      }
      f << str << bc_name.left(10) << "\n";
      bndid++;
    }
  bcid++;
  }
}

void DolfynWriter::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".vrt");
    if (isValid()) {
      QFileInfo filename(getFileName());
      m_VrtFile = new QFile(getFileName());
      m_VrtFile->open(QIODevice::WriteOnly | QIODevice::Text);
      m_CelFile = new QFile(filename.dir().filePath(filename.completeBaseName() + ".cel"));
      m_CelFile->open(QIODevice::WriteOnly | QIODevice::Text);
      m_BndFile = new QFile(filename.dir().filePath(filename.completeBaseName() + ".bnd"));
      m_BndFile->open(QIODevice::WriteOnly | QIODevice::Text);
      
      writeVertices();
      writeElements();
      writeBoundaries();
    
      delete m_VrtFile;
      delete m_CelFile;
      delete m_BndFile;
    }
  } catch (Error err) {
    err.display();
  }
}

