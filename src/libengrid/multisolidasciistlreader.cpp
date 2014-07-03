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

#include "multisolidasciistlreader.h"
#include "stlreader.h"
#include "guimainwindow.h"

#include <QFileInfo>
#include <QInputDialog>
#include <QFile>
#include <QTextStream>


MultiSolidAsciiStlReader::MultiSolidAsciiStlReader()
{
  setFormat("STL files(*.stl *.STL)");
}

void MultiSolidAsciiStlReader::operate()
{
  QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
  readInputFileName(file_info.completeBaseName() + ".stl");
  if (isValid()) {
    double tol = QInputDialog::getText(NULL, "enter STL tolerance", "tolerance", QLineEdit::Normal, "1e-10").toDouble();
    QList<QString> buffer;
    {
      QFile file(getFileName());
      if (!file.open(QFile::ReadOnly)) {
        EG_ERR_RETURN("unable to open file");
      }
      QTextStream f(&file);
      QString buf = "";
      while (!f.atEnd()) {
        QString line = f.readLine();
        buf += line; // endline??
        if (line.left(8) == "endsolid") {
          buffer.append(buf);
          buf = "";
        }
      }
    }
    bool first = true;
    int last_bc = 1;
    foreach (QString buf, buffer) {
      QString file_name = "_" + getFileName();
      {
        QFile file(file_name);
        if (!file.open(QFile::WriteOnly)) {
          EG_ERR_RETURN("unable to open file");
        }
        QTextStream f(&file);
        f << buf << endl;
      }
      StlReader stl;
      stl.setTolerance(tol);
      stl.setFileName(file_name);
      EG_VTKSP(vtkUnstructuredGrid, grid);
      stl.setGrid(grid);
      stl();

      // @todo set boundary names
      EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
      for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
        bc->SetValue(id_cell, last_bc);
      }
      ++last_bc;

      if (first) {
        first = false;
        makeCopy(grid, m_Grid);
      } else {
        MeshPartition part1(m_Grid, true);
        MeshPartition part2(grid, true);
        part1.addPartition(part2);
      }
    }
  }
}
