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
#include "nodelayeriterator.h"

void NodeLayerIterator::operate()
{
  getCells();
  getSurfaceCells(boundary_codes, surface_cells, grid);
  getNodesFromCells(surface_cells, surf_nodes, grid);
  createNodeMapping(surf_nodes, _surf_nodes, grid);
  mark2.fill(false, nodes.size());
  
  // create pairs for first layer
  pair.resize(surf_nodes.size());
  for (int i = 0; i < surf_nodes.size(); ++i) {
    pair[i].item1 = -1;
    pair[i].item2 = surf_nodes[i];
    pair[i].terminate = false;
  };
  
  first_layer = true;
  while (pair.size() > 0) {
    mark1.fill(false, nodes.size());
    pass1();
    for (int i = 0; i < pair.size(); ++i) {
      if (!pair[i].terminate) {
        mark1[_nodes[pair[i].item2]] = true;
      };
    };
    for (int i = 0; i < pair.size(); ++i) {
      if (pair[i].terminate) {
        mark1[_nodes[pair[i].item2]] = false;
      };
    };
    {
      int N = 0;
      for (int i = 0; i < nodes.size(); ++i) {
        if (mark1[i]) {
          ++N;
        };
      };
      item.resize(N);
    };
    {
      int j = 0;
      for (int i = 0; i < nodes.size(); ++i) {
        if (mark1[i]) {
          item[j] = nodes[i];
          ++j;
        };
      };
    };
    pass2();
    foreach (vtkIdType nodeId, item) {
      mark2[_nodes[nodeId]] = true;
    };
    {
      int N = 0;
      for (int i = 0; i < item.size(); ++i) {
        foreach (int j, n2n[_nodes[item[i]]]) {
          if (!mark2[j]) {
            ++N;
          };
        };
      };
      pair.resize(N);
      N = 0;
      for (int i = 0; i < item.size(); ++i) {
        foreach (int j, n2n[_nodes[item[i]]]) {
          if (!mark2[j]) {
            pair[N].item1 = item[i];
            pair[N].item2 = nodes[j];
            pair[N].terminate = false;
            ++N;
          };
        };
      };
    };
    first_layer = false;
  };
};
