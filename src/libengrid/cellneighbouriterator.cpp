// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include "cellneighbouriterator.h"

CellNeighbourIterator::CellNeighbourIterator() 
{
  start = 0; 
}


void CellNeighbourIterator::operate()
{
  getCells();
  pair.resize(1);
  pair[0].item1 = -1;
  pair[0].item2 = start;
  if (_cells[pair[0].item2] < 0) {
    pair[0].item2 = cells[0];
  }
  pair[0].terminate = false;
  mark2.fill(false, cells.size());
  bool first = true;
  while (pair.size() > 0) {
    mark1.fill(false, cells.size());
    if (!first) {
      pass1();
    } else {
      first = false;
    }
    for (int i = 0; i < pair.size(); ++i) {
      if (!pair[i].terminate) {
        mark1[_cells[pair[i].item2]] = true;
      }
    }
    {
      int N = 0;
      for (int i = 0; i < cells.size(); ++i) {
        if (mark1[i]) {
          ++N;
        }
      }
      item.resize(N);
    }
    {
      int j = 0;
      for (int i = 0; i < cells.size(); ++i) {
        if (mark1[i]) {
          item[j] = cells[i];
          ++j;
        }
      }
    }
    pass2();
    for (int i = 0; i < item.size(); ++i) {
      mark2[_cells[item[i]]] = true;
    }
    {
      int N = 0;
      for (int i = 0; i < item.size(); ++i) {
        for (int j = 0; j < c2c[_cells[item[i]]].size(); ++j) {
          if (c2c[_cells[item[i]]][j] >= 0) {
            if (!mark2[c2c[_cells[item[i]]][j]]) {
              ++N;
            }
          }
        }
      }
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
            }
          }
        }
      }
    }
  }
}



