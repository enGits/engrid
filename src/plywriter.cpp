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
#include "plywriter.h"

#include <vtkPLYWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>

#include <QFileInfo>
#include "guimainwindow.h"

PlyWriter::PlyWriter()
{
  setFormat("Stanford PLY files (*.ply *.PLY)");
  m_AsciiFileType = true;
};

void PlyWriter::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".ply");
    if (isValid()) {
      EG_VTKSP(vtkGeometryFilter, geometry);
      geometry->SetInput(m_Grid);
      EG_VTKSP(vtkTriangleFilter, triangle);
      triangle->SetInput(geometry->GetOutput());
      
      EG_VTKSP(vtkPLYWriter, write_ply);
      write_ply->SetInput(triangle->GetOutput());
      write_ply->SetFileName(qPrintable(getFileName()));
      if(m_AsciiFileType) {
        write_ply->SetFileTypeToASCII();
      }
      else {
        write_ply->SetFileTypeToBinary();
      }
      write_ply->SetDataByteOrderToLittleEndian();
      write_ply->SetColorModeToUniformCellColor();
      write_ply->SetColor(255,0,0);
      write_ply->Write();
      
    };
  } catch (Error err) {
    err.display();
  };
};
