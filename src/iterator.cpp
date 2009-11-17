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
#include "iterator.h"

Iterator::Iterator() 
{
  custom_iteration = false;
  volume_iteration = false;
};

void Iterator::pass1()
{
  for (int i = 0; i < pair.size(); ++i) {
    pair[i].terminate = false;
  };
};

void Iterator::getCells()
{
  if (!custom_iteration) {
    if (volume_iteration) {
      getAllVolumeCells(cells, m_Grid);
    } else {
      getAllSurfaceCells(cells, m_Grid);
    };
  };
  getNodesFromCells(cells, nodes, m_Grid);
  createCellMapping(cells, _cells, m_Grid);
  createNodeMapping(nodes, _nodes, m_Grid);
  createNodeToCell(cells, nodes, _nodes, n2c, m_Grid);
  createNodeToNode(cells, nodes, _nodes, n2n, m_Grid);
  createCellToCell(cells, c2c, m_Grid);
};
