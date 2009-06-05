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

#ifndef MESHPARTITION_H
#define MESHPARTITION_H

class MeshPartition;

#include "egvtkobject.h"

class MeshPartition : public EgVtkObject
{

private: // attributes

  /// the grid underlying this mesh partition
  vtkUnstructuredGrid *grid;

  /// all cells of the mesh partition
  QVector<vtkIdType> cells;

  /// inverse indexing for the cells
  QVector<int> _cells;

  /// all nodes of the mesh partition
  QVector<vtkIdType> nodes;

  /// inverse indexing for the nodes
  QVector<int> _nodes;

  /// node to cell information
  QVector<QVector<int> > n2c;

  /// node to node information
  QVector<QVector<int> > n2n;

  /// cell to cell information
  QVector<QVector<int> > c2c;

private: // methods

  /// update x2x structures
  void updateStructures();

public: // methods

  /// Create an empty (undefined) mesh partition
  MeshPartition();

  /**
   * Create a mesh partition with the grid set. Optionally all cells can be selected.
   * @param grid the grid to use
   * @param use_all_cells if set to true all cells will be selected;
   */
  MeshPartition(vtkUnstructuredGrid *grid, bool use_all_cells = false);

  /**
   * Create a mesh partition from a global volume definition
   * @param volume_name the name of the volume
   */
  MeshPartition(QString volume_name);

  /**
   * Set the grid.
   * @param a pointer to the grid
   */
  void setGrid(vtkUnstructuredGrid *grid) { this->grid = grid; }

  /**
   * Access to the grid.
   * @return a pointer to the grid
   */
  vtkUnstructuredGrid* getGrid() const { return grid; }

  /**
   * Define the mesh partition by defining its cells.
   * @param cls the cells of the subset
   */
  template <class C>
  void setCells(const C& cls);

  /**
   * Define the mesh partition by giving a symbolic volume name.
   * The grid will be changed to the default (main) grid that is currently loaded into ENGRID.
   * @param volume_name the symbolic volume name
   */
  void setVolume(QString volume_name);

  /**
   * Define the mesh partition as the remainder of an existing partition.
   * @param part the existing partition
   */
  void setRemainder(const MeshPartition& part);

  /// Access to the cell indices;
  const QVector<vtkIdType>& getCells() const { return cells; }

  /// Access to the local cell indices;
  const QVector<int>& getLocalCells() const { return _cells; }

  /// Access to the node indices;
  const QVector<vtkIdType>& getNodes() const { return nodes; }

  /// Access to the local node indices;
  const QVector<int>& getLocalNodes() const { return _nodes; }

  /// Access to the local node to node structure
  const QVector<QVector<int> >& getN2N() const { return n2n; }

  /// Access to the local node to cell structure
  const QVector<QVector<int> >& getN2C() const { return n2c; }

  /// Access to the local cell to cell structure
  const QVector<QVector<int> >& getC2C() const { return c2c; }

  /// change the face orientation to match the volume definition
  void setVolumeOrientation();

  /// change the orientation to match the original orientation
  void setOriginalOrientation();

  /**
   * Copy the partition to a VTK grid.
   * @param new_grid the grid to copy the partition to (will be resized accordingly).
   */
  void extractToVtkGrid(vtkUnstructuredGrid *new_grid);

  /**
   * Add another partition to this one.
   * At the moment overlapping partitions on two fifferent grids will not be handled well.
   * If both partitions do not have the same underlying grid the grid will be extended in order
   * to add the other partition.
   * @param part the partition to add
   */
  void addPartition(const MeshPartition& part);

  /**
   * compute the smallest edge length of the partition
   * @return the smalles edge length
   */
  double getSmallestEdgeLength() const;

  int n2nLSize(int i_nodes) { return n2n[i_nodes].size(); }
  int n2nLL(int i_nodes, int j) { return n2n[i_nodes][j]; }
  vtkIdType n2nLG(int i_nodes, int j) { return nodes[n2n[i_nodes][j]]; }

  int n2nGSize(vtkIdType id_node) { return n2n[_nodes[id_node]].size(); }
  int n2nGL(vtkIdType id_node, int j) { return n2n[_nodes[id_node]][j]; }
  vtkIdType n2nGG(vtkIdType id_node, int j) { return nodes[n2n[_nodes[id_node]][j]]; }

  int n2cLSize(int i_nodes) { return n2c[i_nodes].size(); }
  int n2cLL(int i_nodes, int j) { return n2c[i_nodes][j]; }
  vtkIdType n2cLG(int i_nodes, int j) { return nodes[n2c[i_nodes][j]]; }

  int n2cGSize(vtkIdType id_node) { return n2c[_nodes[id_node]].size(); }
  int n2cGL(vtkIdType id_node, int j) { return n2n[_nodes[id_node]][j]; }
  vtkIdType n2cGG(vtkIdType id_node, int j) { return cells[n2n[_nodes[id_node]][j]]; }

  int c2cLSize(int i_cells) { return c2c[i_cells].size(); }
  int c2cLL(int i_cells, int j) { return c2c[i_cells][j]; }
  vtkIdType c2cLG(int i_cells, int j) { return cells[c2c[i_cells][j]]; }

  int c2cGSize(vtkIdType id_cell) { return c2c[_cells[id_cell]].size(); }
  int c2cGL(vtkIdType id_cell, int j) { return c2c[_cells[id_cell]][j]; }
  vtkIdType c2cGG(vtkIdType id_cell, int j) { return cells[c2c[_cells[id_cell]][j]]; }

};



template <class C>
void MeshPartition::setCells(const C& cls)
{
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  updateStructures();
}

#endif // MESHPARTITION_H
