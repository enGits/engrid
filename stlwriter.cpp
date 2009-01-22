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
#include "stlwriter.h"

#include <vtkSTLWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>

StlWriter::StlWriter()
{
  setFormat("stereolithography files (*.stl *.STL)");
};

void StlWriter::operate()
{
  try {
    readOutputFileName();
    if (isValid()) {
      EG_VTKSP(vtkGeometryFilter, geometry);
      geometry->SetInput(grid);
      EG_VTKSP(vtkTriangleFilter, triangle);
      triangle->SetInput(geometry->GetOutput());
      EG_VTKSP(vtkSTLWriter, write_stl);
      write_stl->SetInput(triangle->GetOutput());
      write_stl->SetFileName(getFileName().toAscii().data());
      write_stl->SetFileTypeToASCII();
      write_stl->Write();
    };
  } catch (Error err) {
    err.display();
  };
};


