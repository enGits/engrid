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
#include "vtkEgExtractVolumeCells.h"

vtkStandardNewMacro(vtkEgExtractVolumeCells);

vtkEgExtractVolumeCells::vtkEgExtractVolumeCells()
{
  SetClippingOff();
  SetX(vec3_t(0,0,0));
  SetN(vec3_t(0,0,1));
  tetra   = true;
  pyramid = true;
  wedge   = true;
  hexa    = true;
};

void vtkEgExtractVolumeCells::ExecuteEg()
{
  QSet<vtkIdType> ex_cells;
  for (vtkIdType id_cell = 0; id_cell < input->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, input)) {
      bool select = true;
      vtkIdType type_cell = input->GetCellType(id_cell);
      if (!tetra   && type_cell == VTK_TETRA)      select = false;
      if (!pyramid && type_cell == VTK_PYRAMID)    select = false;
      if (!wedge   && type_cell == VTK_WEDGE)      select = false;
      if (!hexa    && type_cell == VTK_HEXAHEDRON) select = false;
      if (Clip && select) {
        vtkIdType *pts;
        vtkIdType  N_pts;
        input->GetCellPoints(id_cell, N_pts, pts);
        for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
          vec3_t x;
          input->GetPoints()->GetPoint(pts[i_pts],x.data());
          if ((x-X)*N < 0) {
            select = false;
            break;
          };
        };
      };
      if (select) {
        ex_cells.insert(id_cell);
      };
    };
  };
  QVector<vtkIdType> cells(ex_cells.size());
  qCopy(ex_cells.begin(), ex_cells.end(), cells.begin());
  QVector<vtkIdType> nodes;
  QVector<int>       _nodes;
  getNodesFromCells(cells, nodes, input);
  createNodeMapping(nodes, _nodes, input);
  allocateGrid(output, cells.size(), nodes.size());
  vtkIdType cellId, nodeId;
  foreach(nodeId, nodes) {
    vec3_t x;
    input->GetPoints()->GetPoint(nodeId, x.data());
    output->GetPoints()->SetPoint(_nodes[nodeId], x.data());
  };
  foreach(cellId, cells) {
    vtkIdType  Npts;
    vtkIdType *pts;
    input->GetCellPoints(cellId, Npts, pts);
    vtkIdType *new_pts = new vtkIdType[Npts];
    for (int i = 0; i < Npts; ++i) {
      new_pts[i] = _nodes[pts[i]];
    };
    vtkIdType newId = output->InsertNextCell(input->GetCellType(cellId), Npts, new_pts);
    copyCellData(input, cellId, output, newId);
    delete [] new_pts;
  };
};

