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
#include "stlwriter.h"

#include <vtkSTLWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>

#include <QFileInfo>
#include "guimainwindow.h"

StlWriter::StlWriter()
{
  setFormat("stereolithography files (*.stl *.STL)");
  m_AsciiFileType = true;
};

void StlWriter::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readOutputFileName(file_info.completeBaseName() + ".stl");
    if (isValid()) {
      EG_VTKSP(vtkGeometryFilter, geometry);
      geometry->SetInputData(m_Grid);
      EG_VTKSP(vtkTriangleFilter, triangle);
      triangle->SetInputConnection(geometry->GetOutputPort());
      EG_VTKSP(vtkSTLWriter, write_stl);
      write_stl->SetInputConnection(  triangle->GetOutputPort());
      write_stl->SetFileName(qPrintable(getFileName()));
      if(m_AsciiFileType) {
        write_stl->SetFileTypeToASCII();
      }
      else {
        write_stl->SetFileTypeToBinary();
      }
      write_stl->Write();
    };
  } catch (Error err) {
    err.display();
  };
};
