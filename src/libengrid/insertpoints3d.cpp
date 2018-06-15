// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2018 enGits GmbH                                      +
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
#include "insertpoints3d.h"

#include <vtkCharArray.h>

void InsertPoints3d::reset()
{
  m_EdgeBuffer.clear();
}

bool InsertPoints3d::markEdge(vtkIdType id_node1, vtkIdType id_node2, int type, vec3_t x)
{
  stencil_t E;
  E.p1 = id_node1;
  E.p2 = id_node2;
  getEdgeCells(id_node1, id_node2, E.id_cell);
  QSet<vtkIdType> nodes;
  foreach (vtkIdType id_cell, E.id_cell) {
    if (m_CellMarked[m_Part.localCell(id_cell)]) {
      return false;
    }
  }
  foreach (vtkIdType id_cell, E.id_cell) {
    EG_GET_CELL_AND_TYPE(id_cell, m_Grid);
    if (type_cell != VTK_TETRA) {
      EG_BUG;
    }
    m_CellMarked[m_Part.localCell(id_cell)] = true;
    for (vtkIdType i = 0; i < num_pts; ++i) {
      nodes.insert(pts[i]);
    }
  }
  E.id_node = set2Vector(nodes, false);
  m_Edges << E;
  m_X     << x;
  m_Type  << type;
  return true;
}

bool InsertPoints3d::addEdge(vtkIdType id_node1, vtkIdType id_node2, int type, vec3_t x)
{
  edge_t E(id_node1, id_node2, type, x);
  if (m_EdgeBuffer.contains(E)) {
    return false;
  }
  m_EdgeBuffer << E;
  return true;
}

bool InsertPoints3d::addEdge(vtkIdType id_node1, vtkIdType id_node2, int type)
{
  vec3_t x1, x2;
  m_Grid->GetPoint(id_node1, x1.data());
  m_Grid->GetPoint(id_node2, x2.data());
  return addEdge(id_node1, id_node2, type, 0.5*(x1 + x2));
}

void InsertPoints3d::getOrderedCellNodes(stencil_t E, vtkIdType id_cell, QVector<vtkIdType> &ordered_nodes)
{
  EG_GET_CELL(id_cell, m_Grid);
  QList<vtkIdType> nodes;
  for (vtkIdType i = 0; i < num_pts; ++i) {
    nodes << pts[i];
  }
  nodes << nodes;
  bool started = false;
  ordered_nodes.clear();
  foreach (vtkIdType id_node, nodes) {
    if (id_node == E.p1 || id_node == E.p2) {
      started = true;
    }
    if (started) {
      if (ordered_nodes.size() < num_pts) {
        ordered_nodes << id_node;
      } else {
        break;
      }
    }
  }
  if (ordered_nodes.size() != num_pts) {
    EG_BUG;
  }
}

void InsertPoints3d::getSplitCellNodes(stencil_t E, vtkIdType id_cell, vtkIdType id_new_node, QVector<vtkIdType> &nodes1, QVector<vtkIdType> &nodes2)
{
  // this only works for tetras!!
  //
  QVector<vtkIdType> ordered_nodes;
  getOrderedCellNodes(E, id_cell, ordered_nodes);
  nodes1.clear();
  nodes2.clear();
  //
  if        (ordered_nodes[2] == E.p1 || ordered_nodes[2] == E.p2) {
    //
    // case 1
    //
    nodes1 << ordered_nodes[0] << id_new_node << ordered_nodes[3] << ordered_nodes[1];
    nodes2 << id_new_node << ordered_nodes[2] << ordered_nodes[3] << ordered_nodes[1];
    //
  } else if (ordered_nodes[1] == E.p1 || ordered_nodes[1] == E.p2) {
    //
    // case 2
    //
    nodes1 << ordered_nodes[0] << id_new_node << ordered_nodes[2] << ordered_nodes[3];
    nodes2 << id_new_node << ordered_nodes[1] << ordered_nodes[2] << ordered_nodes[3];
    //
  } else if (ordered_nodes[3] == E.p1 || ordered_nodes[3] == E.p2) {
    //
    // case 3
    //
    nodes1 << ordered_nodes[0] << ordered_nodes[2] << id_new_node << ordered_nodes[1];
    nodes2 << id_new_node << ordered_nodes[2] << ordered_nodes[3] << ordered_nodes[1];
    //
  } else {
    EG_BUG;
  }
}

bool InsertPoints3d::cellMarked(vtkIdType id_cell)
{
  int i_cell = m_Part.localCell(id_cell);
  if (i_cell < 0) {
    return false;
  }
  return m_CellMarked[i_cell];
}

void InsertPoints3d::splitIteration()
{
  // determine new number of nodes and cells and allocate new grid
  //
  int num_new_cells = m_Grid->GetNumberOfCells();
  int num_new_nodes = m_Grid->GetNumberOfPoints();
  foreach (stencil_t E, m_Edges) {
    num_new_nodes += 1;
    num_new_cells += E.id_cell.size();
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, num_new_cells, num_new_nodes);
  //
  // copy all unaffected cells
  //
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (!cellMarked(id_cell)) {
      vtkIdType id_new_cell = copyCell(m_Grid, id_cell, new_grid);
      copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
    }
  }
  //
  // copy all existing nodes
  //
  EG_FORALL_NODES(id_node, m_Grid) {
    vec3_t x;
    m_Grid->GetPoints()->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(m_Grid, id_node, new_grid, id_node);
  }
  //
  // create new nodes and split cells
  //
  EG_VTKDCN(vtkCharArray, node_type, new_grid, "node_type");
  vtkIdType id_node = m_Grid->GetNumberOfPoints();
  for (int i = 0; i < m_Edges.size(); ++i) {
    new_grid->GetPoints()->SetPoint(id_node, m_X[i].data());
    copyNodeData(m_Grid, m_Edges[i].p1, new_grid, id_node);
    node_type->SetValue(id_node, m_Type[i]);
    QVector<vtkIdType> nodes1, nodes2;
    foreach (vtkIdType id_cell, m_Edges[i].id_cell) {
      //
      // this only works for tetras!
      //
      getSplitCellNodes(m_Edges[i], id_cell, id_node, nodes1, nodes2);
      //
      vtkIdType new_pts1[4], new_pts2[4];
      for (int j = 0; j < 4; ++j) {
        new_pts1[j] = nodes1[j];
        new_pts2[j] = nodes2[j];
      }
      vtkIdType id_cell1 = new_grid->InsertNextCell(VTK_TETRA, 4, new_pts1);
      vtkIdType id_cell2 = new_grid->InsertNextCell(VTK_TETRA, 4, new_pts2);
      copyCellData(m_Grid, id_cell, new_grid, id_cell1);
      copyCellData(m_Grid, id_cell, new_grid, id_cell2);
    }
    ++id_node;
  }
  //
  makeCopy(new_grid, m_Grid);
}

void InsertPoints3d::operate()
{
  while (m_EdgeBuffer.size() > 0) {
    //
    setAllVolumeCells();
    m_Edges.clear();
    m_X.clear();
    m_Type.clear();
    m_CellMarked.fill(false, m_Part.getNumberOfCells());
    //
    QList<edge_t> rest_edges;
    foreach (edge_t E, m_EdgeBuffer) {
      if (!markEdge(E.p1, E.p2, E.type, E.x)) {
        rest_edges << E;
      }
    }
    splitIteration();
    m_EdgeBuffer = rest_edges;
    //m_EdgeBuffer.clear();
  }
}
