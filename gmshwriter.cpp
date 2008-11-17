//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "gmshwriter.h"

GmshWriter::GmshWriter()
{
  setFormat("Gmsh files(*.msh)");
};

void GmshWriter::writeAscii1(vtkUnstructuredGrid *grid)
{
  QFile file(getFileName());
  file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&file);
  f << "$NOD\n";
  f << grid->GetNumberOfPoints() << '\n';
  for (vtkIdType nodeId = 0; nodeId < grid->GetNumberOfPoints(); ++nodeId) {
    vec3_t x;
    grid->GetPoints()->GetPoint(nodeId, x.data());
    f << nodeId+1 << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
  };
  f << "$ENDNOD\n";
  f << "$ELM\n";
  f << grid->GetNumberOfCells() << '\n';
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtkIdType  Npts;
    vtkIdType *pts;
    grid->GetCellPoints(cellId, Npts, pts);
    f << cellId+1;
    if (grid->GetCellType(cellId) == VTK_TRIANGLE) {
      f << " 2 " << cell_code->GetValue(cellId) << " 0 " << Npts;
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_QUAD) {
      f << " 3 " << cell_code->GetValue(cellId) << " 0 " << Npts;
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_TETRA) {
      f << " 4 0 0 " << Npts;
      
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
      
      //f << ' ' << pts[0]+1 << ' ' << pts[1]+1 << ' ' << pts[3]+1 << ' ' << pts[2]+1;
    } else if (grid->GetCellType(cellId) == VTK_PYRAMID) {
      f << " 7 0 0 " << Npts;
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_WEDGE) {
      f << " 6 0 0 " << Npts;
      for (int i = 3; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
      for (int i = 0; i < Npts-3; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_HEXAHEDRON) {
      f << " 5 0 0 " << Npts;
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    };
    f << '\n';
  };
  f << "$ENDELM\n";
};

void GmshWriter::writeAscii2(vtkUnstructuredGrid *grid)
{
  QFile file(getFileName());
  file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&file);
  f << "$MeshFormat\n2.0 0 8\n$EndMeshFormat\n";
  f << "$Nodes\n";
  f << grid->GetNumberOfPoints() << '\n';
  for (vtkIdType nodeId = 0; nodeId < grid->GetNumberOfPoints(); ++nodeId) {
    vec3_t x;
    grid->GetPoints()->GetPoint(nodeId, x.data());
    f << nodeId+1 << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
  };
  f << "$EndNodes\n";
  f << "$Elements\n";
  f << grid->GetNumberOfCells() << '\n';
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
    vtkIdType  Npts;
    vtkIdType *pts;
    grid->GetCellPoints(cellId, Npts, pts);
    f << cellId+1;
    if (grid->GetCellType(cellId) == VTK_TRIANGLE) {
      f << " 2 1 " << cell_code->GetValue(cellId);
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_QUAD) {
      f << " 3 1 " << cell_code->GetValue(cellId);
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_TETRA) {
      f << " 4 1 0 ";
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_PYRAMID) {
      f << " 7 1 0 ";
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_WEDGE) {
      f << " 6 1 0 ";
      for (int i = 3; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
      for (int i = 0; i < Npts-3; ++i) {
        f << ' ' << pts[i]+1;
      };
    } else if (grid->GetCellType(cellId) == VTK_HEXAHEDRON) {
      f << " 5 1 0 ";
      for (int i = 0; i < Npts; ++i) {
        f << ' ' << pts[i]+1;
      };
    };
    f << '\n';
  };
  f << "$EndElements\n";
};

void GmshWriter::operate()
{
  try {
    readOutputFileName();
    if (isValid()) {
      if      (format == ascii1) writeAscii1(grid);
      else if (format == ascii2) writeAscii2(grid);
    };
  } catch (Error err) {
    err.display();
  };
};
