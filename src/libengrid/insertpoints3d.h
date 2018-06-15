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
#ifndef INSERTPOINTS3D_H
#define INSERTPOINTS3D_H

#include "surfaceoperation.h"

class InsertPoints3d;


class InsertPoints3d : public SurfaceOperation
{

protected: // data types

  struct edge_t
  {
    vtkIdType p1, p2;
    int       type;
    vec3_t    x;

    edge_t(vtkIdType a_p1, vtkIdType a_p2, int a_type, vec3_t an_x)
    {
      if (a_p1 > a_p2) {
        p1 = a_p2;
        p2 = a_p1;
      } else {
        p1 = a_p1;
        p2 = a_p2;
      }
      type = a_type;
      x = an_x;
    }

    bool operator== (const edge_t E)
    {
      if (p1 == E.p1 && p2 == E.p2) return true;
      return false;
    }
  };


protected: // attributes

  QList<edge_t>    m_EdgeBuffer;
  QList<stencil_t> m_Edges;
  QList<vec3_t>    m_X;
  QList<int>       m_Type;
  QVector<bool>    m_CellMarked;


protected: // methods

  virtual void operate();

  void getOrderedCellNodes(stencil_t E, vtkIdType id_cell, QVector<vtkIdType> &ordered_nodes);
  void getSplitCellNodes  (stencil_t E, vtkIdType id_cell, vtkIdType id_new_node, QVector<vtkIdType> &nodes1, QVector<vtkIdType> &nodes2);

  bool markEdge(vtkIdType id_node1, vtkIdType id_node2, int type, vec3_t x);
  bool cellMarked(vtkIdType id_cell);
  void splitIteration();


public: // methods

  void reset();

  /**
   * @brief Add and edge for splitting.
   * @param id_node1 first node of the edge
   * @param id_node2 second node of the edge
   * @param type the type of the new node
   * @param x coordinates of the new node
   * @return true if the edge has not been added before
   */
  bool addEdge(vtkIdType id_node1, vtkIdType id_node2, int type, vec3_t x);

  /**
   * @brief Add and edge for splitting.
   * @param id_node1 first node of the edge
   * @param id_node2 second node of the edge
   * @param type the type of the new node
   * @return true if the edge has not been added before
   */
  bool addEdge(vtkIdType id_node1, vtkIdType id_node2, int type);


};

#endif // INSERTPOINTS3D_H
