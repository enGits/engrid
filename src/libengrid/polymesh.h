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
#ifndef polymesh_H
#define polymesh_H

class PolyMesh;

#include "egvtkobject.h"
#include "meshpartition.h"

class PolyMesh : public EgVtkObject
{
  
protected: // data types

  struct face_t {
    QVector<int> node;
    int owner, neighbour;
    int bc;
    int operator[](int i);
    bool operator<(const face_t &F) const;
    vec3_t normal(vtkUnstructuredGrid *grid);
    face_t() {}
    face_t(int N, int o, int n, int b = 0);
  };

  struct node_t {
    QVector<vtkIdType> id;
    node_t(vtkIdType id1);
    node_t(vtkIdType id1, vtkIdType id2);
    node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3);
    node_t(vtkIdType id1, vtkIdType id2, vtkIdType id3, vtkIdType id4);
    bool operator<(const node_t &N) const;
    bool operator>(const node_t &N) const;
    bool operator==(const node_t &N) const;
  };

  
protected: // attributes

  vtkUnstructuredGrid *m_Grid;
  MeshPartition        m_Part;
  QVector<int>         m_Cell2PCell;
  QVector<int>         m_Node2PCell;
  QVector<node_t>      m_Nodes;
  QVector<face_t>      m_Faces;

protected: // methods

  bool isHexCoreNode(vtkIdType) { return false; }
  bool isHexCoreCell(vtkIdType) { return false; }

  /**
    * Get internal indices of adjacent faces of an edge.
    * See http://engits.eu/wiki/index.php/Manual/Element_Types for details about the internal indices.
    * @param id_cell global id of the cell (only tetras allowed)
    * @param id_node1 global id of the first node of the edge
    * @param id_node2 global id of the second node of the edge
    * @param face1 (return value) will hold the internal index of the first face (0,1,2,3)
    * @param face2 (return value) will hold the internal index of the second face (0,1,2,3)
    */
  void getFacesOfEdgeInsideTetra(vtkIdType id_cell, vtkIdType id_node1, vtkIdType id_node2, int &face1, int &face2);

  /**
    * Find all cells around an edge.
    * This method will find all tetra cells around an edge in the mesh.
    * It can either be a closed loop of tetras or a list which is terminated by surface elements.
    * @param id_node1 first node of the edge
    * @param id_node2 second node of the edge
    * @param cells list of tetras/faces around the edge
    */
  void getSortedEdgeCells(vtkIdType id_node1, vtkIdType id_node2, QList<vtkIdType> &cells);

  void findDualCells();
  void countNodesAndFaces();
   
public: // methods
  
  PolyMesh(vtkUnstructuredGrid *grid, bool dual_mesh = true);
  int totalNumNodes() const         { EG_BUG; return 0; }
  vec3_t nodeVector(int i) const    { EG_BUG; return vec3_t(0,0,0); }
  int numNodes(int i) const         { EG_BUG; return 0; }
  int nodeIndex(int i, int j) const { EG_BUG; return 0; }
  int numFaces() const              { EG_BUG; return 0; }
  //int numCells() const { return id_newpc; }
  int owner(int i) const            { EG_BUG; return 0; }
  int neighbour(int i) const        { EG_BUG; return 0; }
  int boundaryCode(int i) const     { EG_BUG; return 0; }
  int numBCs() const                { EG_BUG; return 0; }
  //void getFace(vtkIdType idx, int subidx, QList<vtkIdType> &nodes);
  //void getEdge(vtkIdType idx, int subidx, QList<vtkIdType> &nodes);
  
  
};

#endif
