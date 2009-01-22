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
#include "neutralwriter.h"

NeutralWriter::NeutralWriter()
{
  setFormat("Neutral mesh files(*.mesh)");
};

void NeutralWriter::operate()
{
  try {
    readOutputFileName();
    if (isValid()) {
      QFile file(getFileName());
      file.open(QIODevice::WriteOnly | QIODevice::Text);
      QTextStream f(&file);
      f << grid->GetNumberOfPoints() << "\n";
      for (vtkIdType pointId = 0; pointId < grid->GetNumberOfPoints(); ++pointId) {
        vec3_t x;
        grid->GetPoints()->GetPoint(pointId, x.data());
        f << x[0] << " " << x[1] << " " << x[2] << "\n";
      };
      vtkIdType Nvol = 0;
      vtkIdType Nsurf = 0;
      for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
        vtkIdType cellType = grid->GetCellType(cellId);
        if ((cellType != VTK_TRIANGLE) && (cellType != VTK_TETRA)) {
          EG_ERR_RETURN("only simplex elements are allowed for the NEUTRAL format");
        };
        if (isSurface(cellId, grid)) {
          ++Nsurf;
        } else {
          ++Nvol;
        };
      };
      f << Nvol << "\n";
      for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
        if (!isSurface(cellId, grid)) {
          vtkIdType Npts, *pts;
          grid->GetCellPoints(cellId, Npts, pts);
          f << "1 " << pts[0]+1 << " " << pts[1]+1 << " " << pts[3]+1 << " " << pts[2]+1 << "\n";
        };
      };
      f << Nsurf << "\n";
      EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
      for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
        if (isSurface(cellId, grid)) {
          vtkIdType Npts, *pts;
          grid->GetCellPoints(cellId, Npts, pts);
          f << 1;//cell_code->GetValue(cellId);
          f << " " << pts[2]+1 << " " << pts[1]+1 << " " << pts[0]+1 << "\n";
        };
      };
    };
  } catch (Error err) {
    err.display();
  };
};
