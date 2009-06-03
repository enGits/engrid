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

#include "egvtkobject.h"

#include "guimainwindow.h"

class MeshPartition : public EgVtkObject
{

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
   * Define the mesh partition as the remainder of an existing partition.
   * @param part the existing partition
   */
  void setRemainder(const MeshPartition& part);

  /// Access to the cell indices;
  const QVector<vtkIdType>& getCells() { return cells; }

  /// Access to the local cell indices;
  const QVector<int>& getLocalCells() { return _cells; }

  /// Access to the node indices;
  const QVector<vtkIdType>& getNodes() { return nodes; }

  /// Access to the local node indices;
  const QVector<int>& getLocalNodes() { return _nodes; }

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

};



template <class C>
void MeshPartition::setCells(const C& cls)
{
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  createCellMapping(cells, _cells, grid);
  createNodeMapping(nodes, _nodes, grid);
}

#endif // MESHPARTITION_H
