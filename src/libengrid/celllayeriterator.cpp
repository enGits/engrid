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
#include "celllayeriterator.h"

void CellLayerIterator::operate()
{
  QVector<vtkIdType> surf_cells, surf_nodes;
  QVector<int>       _surf_nodes;
  
  custom_iteration = true;
  getAllCells(cells, m_Grid);
  getCells();
  
  getSurfaceCells(m_BoundaryCodes, surf_cells, m_Grid);
  getNodesFromCells(surf_cells, surf_nodes, m_Grid);
  createNodeMapping(surf_nodes, _surf_nodes, m_Grid);
  
  // create first layer
  item.resize(surf_cells.size());
  pair.resize(surf_cells.size());
  for (int i_scell = 0; i_scell < surf_cells.size(); ++i_scell) {
    vtkIdType id_cell1 = surf_cells[i_scell];
    pair[i_scell].item1 = id_cell1;
    EG_GETPTS(pts, id_cell1, m_Grid);
    bool ok = false;
    if (Npts > 0) {
      foreach (int i_cells_a, n2c[_nodes[pts[0]]]) {
        vtkIdType id_cell_a = cells[i_cells_a];
        if (!isSurface(id_cell_a, m_Grid)) {
          ok = true;
          for (int i = 1; i < Npts; ++i) {
            if (ok) {
              ok = false;
              foreach (int i_cells_b, n2c[_nodes[pts[i]]]) {
                vtkIdType id_cell_b = cells[i_cells_b];
                if (id_cell_a == id_cell_b) {
                  ok = true;
                  break;
                };
              };
            };
          };
        };
        if (ok) {
          pair[i_scell].item2 = id_cell_a;
          pair[i_scell].terminate = false;
          break;
        };
      };
    };
    if (!ok) {
      EG_BUG;
    };
  };
  
  mark2.fill(false, cells.size());
  
  first_layer = true;
  
  while (pair.size() > 0) {
    mark1.fill(false, cells.size());
    pass1();
    for (int i = 0; i < pair.size(); ++i) {
      if (!pair[i].terminate) {
        mark1[_cells[pair[i].item2]] = true;
      };
    };
    for (int i = 0; i < pair.size(); ++i) {
      if (pair[i].terminate) {
        mark1[_cells[pair[i].item2]] = false;
      };
    };
    {
      int N = 0;
      for (int i = 0; i < cells.size(); ++i) {
        if (mark1[i]) {
          ++N;
        };
      };
      item.resize(N);
    };
    {
      int j = 0;
      for (int i = 0; i < cells.size(); ++i) {
        if (mark1[i]) {
          item[j] = cells[i];
          ++j;
        };
      };
    };
    pass2();
    for (int i = 0; i < item.size(); ++i) {
      mark2[_cells[item[i]]] = true;
    };
    {
      int N = 0;
      for (int i = 0; i < item.size(); ++i) {
        for (int j = 0; j < c2c[_cells[item[i]]].size(); ++j) {
          if (c2c[_cells[item[i]]][j] >= 0) {
            if (!mark2[c2c[_cells[item[i]]][j]]) {
              ++N;
            };
          };
        };
      };
      pair.resize(N);
      N = 0;
      for (int i = 0; i < item.size(); ++i) {
        for (int j = 0; j < c2c[_cells[item[i]]].size(); ++j) {
          if (c2c[_cells[item[i]]][j] >= 0) {
            if (!mark2[c2c[_cells[item[i]]][j]]) {
              pair[N].item1 = item[i];
              pair[N].item2 = cells[c2c[_cells[item[i]]][j]];
              pair[N].terminate = false;
              ++N;
            };
          };
        };
      };
    };
    first_layer = false;
  };
};


