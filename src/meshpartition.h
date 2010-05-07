//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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


  vtkUnstructuredGrid*   m_Grid;   ///< the grid underlying this mesh partition
  QVector<vtkIdType>     m_Cells;  ///< all cells of the mesh partition
  QVector<int>           m_LCells; ///< inverse indexing for the cells
  QVector<vtkIdType>     m_Nodes;  ///< all nodes of the mesh partition
  QVector<int>           m_LNodes; ///< inverse indexing for the nodes
  QVector<QVector<int> > m_N2C;    ///< node to cell information
  QVector<QVector<int> > m_N2N;    ///< node to node information
  QVector<QVector<int> > m_C2C;    ///< cell to cell information

  int m_CellsStamp;  ///< "time"-stamp
  int m_LCellsStamp; ///< "time"-stamp
  int m_NodesStamp;  ///< "time"-stamp
  int m_LNodesStamp; ///< "time"-stamp
  int m_N2NStamp;    ///< "time"-stamp
  int m_N2CStamp;    ///< "time"-stamp
  int m_C2CStamp;    ///< "time"-stamp

private: // methods

  void resetTimeStamps();
  void checkNodes();
  void checkLCells();
  void checkLNodes();
  void checkN2N();
  void checkN2C();
  void checkC2C();

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
  void setGrid(vtkUnstructuredGrid *grid) { m_Grid = grid; }

  /**
   * Access to the grid.
   * @return a pointer to the grid
   */
  vtkUnstructuredGrid* getGrid() const { return m_Grid; }

  /**
   * Define the mesh partition by defining its cells.
   * @param cls the cells of the subset
   */
  template <class C>
  void setCells(const C& cls);

  /**
   * Define the mesh partition by defining all its cells.
   */
  void setAllCells();
  
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

  const QVector<vtkIdType>&     getCells() const; ///< Access to the cell indices
  const QVector<int>&           getLocalCells();  ///< Access to the local cell indices
  const QVector<vtkIdType>&     getNodes();       ///< Access to the node indices
  const QVector<int>&           getLocalNodes();  ///< Access to the local node indices
  const QVector<QVector<int> >& getN2N();         ///< Access to the local node to node structure
  const QVector<QVector<int> >& getN2C();         ///< Access to the local node to cell structure
  const QVector<QVector<int> >& getC2C();         ///< Access to the local cell to cell structure

  void setVolumeOrientation();   ///< change the face orientation to match the volume definition
  void setOriginalOrientation(); ///< change the orientation to match the original orientation

  /**
   * Copy the partition to a VTK grid.
   * @param new_grid the grid to copy the partition to (will be resized accordingly).
   */
  void extractToVtkGrid(vtkUnstructuredGrid *new_grid);

  /**
   * Add another partition to this one.
   * At the moment overlapping partitions on two different grids will not be handled well.
   * If both partitions do not have the same underlying grid the grid will be extended in order
   * to add the other partition.
   * @param part the partition to add
   */
  void addPartition(const MeshPartition& part);

  /**
   * compute the smallest edge length of the partition
   * @return the smallest edge length
   */
  double getSmallestEdgeLength() const;

  /**
   * Get the number of nodes in the partition.
   * @return the number of nodes
   */
  int getNumberOfNodes();

  /**
   * Get the number of cells in the partition.
   * @return the number of cells
   */
  int getNumberOfCells();

  int       localNode(vtkIdType id_node);
  vtkIdType globalNode(int i);
  int       localCell(vtkIdType id_cell);
  vtkIdType globalCell(int i);


  int       n2nLSize(int i_nodes);
  int       n2nLL(int i_nodes, int j);
  vtkIdType n2nLG(int i_nodes, int j);
  int       n2nGSize(vtkIdType id_node);
  int       n2nGL(vtkIdType id_node, int j);
  vtkIdType n2nGG(vtkIdType id_node, int j);
  int       n2cLSize(int i_nodes);
  int       n2cLL(int i_nodes, int j);
  vtkIdType n2cLG(int i_nodes, int j);
  int       n2cGSize(vtkIdType id_node);
  int       n2cGL(vtkIdType id_node, int j);
  vtkIdType n2cGG(vtkIdType id_node, int j);
  int       c2cLSize(int i_cells);
  int       c2cLL(int i_cells, int j);
  vtkIdType c2cLG(int i_cells, int j);
  int       c2cGSize(vtkIdType id_cell);
  int       c2cGL(vtkIdType id_cell, int j);
  vtkIdType c2cGG(vtkIdType id_cell, int j);

};


template <class C>
inline void MeshPartition::setCells(const C& cls)
{
  m_Cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), m_Cells.begin());
  ++m_CellsStamp;
}

inline void MeshPartition::setAllCells() {
  QVector <vtkIdType> all_cells;
  getAllCells(all_cells, m_Grid);
  this->setCells(all_cells);
}

inline void MeshPartition::checkNodes()
{
  if (m_CellsStamp > m_NodesStamp) {
    getNodesFromCells(m_Cells, m_Nodes, m_Grid);
    m_NodesStamp = m_CellsStamp;
  }
}

inline void MeshPartition::checkLCells()
{
  if (m_CellsStamp > m_LCellsStamp) {
    createCellMapping(m_Cells, m_LCells, m_Grid);
    m_LCellsStamp = m_CellsStamp;
  }
}

inline void MeshPartition::checkLNodes()
{
  checkNodes();
  if (m_NodesStamp > m_LNodesStamp) {
    createNodeMapping(m_Nodes, m_LNodes, m_Grid);
    m_LNodesStamp = m_NodesStamp;
  }
}

inline void MeshPartition::checkN2N()
{
  checkLNodes();
  if (m_LNodesStamp > m_N2NStamp) {
    createNodeToNode(m_Cells, m_Nodes, m_LNodes, m_N2N, m_Grid);
    m_N2NStamp = m_LNodesStamp;
  }
}

inline void MeshPartition::checkN2C()
{
  checkLNodes();
  if (m_LNodesStamp > m_N2CStamp) {
    createNodeToCell(m_Cells, m_Nodes, m_LNodes, m_N2C, m_Grid);
    m_N2CStamp = m_LNodesStamp;
  }
}

inline void MeshPartition::checkC2C()
{
  if (m_CellsStamp > m_C2CStamp) {
    createCellToCell(m_Cells, m_C2C, m_Grid);
    m_C2CStamp = m_CellsStamp;
  }
}

inline const QVector<vtkIdType>& MeshPartition::getCells() const
{
  return m_Cells;
}

inline const QVector<int>& MeshPartition::getLocalCells()
{
  checkLCells();
  return m_LCells;
}

inline const QVector<vtkIdType>& MeshPartition::getNodes()
{
  checkNodes();
  return m_Nodes;
}

inline const QVector<int>& MeshPartition::getLocalNodes()
{
  checkLNodes();
  return m_LNodes;
}

inline const QVector<QVector<int> >& MeshPartition::getN2N()
{
  checkN2N();
  return m_N2N;
}

inline const QVector<QVector<int> >& MeshPartition::getN2C()
{
  checkN2C();
  return m_N2C;
}

inline const QVector<QVector<int> >& MeshPartition::getC2C()
{
  checkC2C();
  return m_C2C;
}

inline int MeshPartition::n2nLSize(int i_nodes)
{
  checkN2N();
  return m_N2N[i_nodes].size();
}

inline int MeshPartition::n2nLL(int i_nodes, int j)
{
  checkN2N();
  return m_N2N[i_nodes][j];
}

inline vtkIdType MeshPartition::n2nLG(int i_nodes, int j)
{
  checkN2N();
  return m_Nodes[m_N2N[i_nodes][j]];
}

inline int MeshPartition::n2nGSize(vtkIdType id_node)
{
  checkN2N();
  return m_N2N[m_LNodes[id_node]].size();
}

inline int MeshPartition::n2nGL(vtkIdType id_node, int j)
{
  checkN2N();
  return m_N2N[m_LNodes[id_node]][j];
}

inline vtkIdType MeshPartition::n2nGG(vtkIdType id_node, int j)
{
  checkN2N();
  return m_Nodes[m_N2N[m_LNodes[id_node]][j]];
}

inline int MeshPartition::n2cLSize(int i_nodes)
{
  checkN2C();
  return m_N2C[i_nodes].size();
}

inline int MeshPartition::n2cLL(int i_nodes, int j)
{
  checkN2C();
  return m_N2C[i_nodes][j];
}

inline vtkIdType MeshPartition::n2cLG(int i_nodes, int j)
{
  checkN2C();
  int i_cell = m_N2C[i_nodes][j];
  if(i_cell<0) return(-1);
  else return m_Cells[i_cell];
}

inline int MeshPartition::n2cGSize(vtkIdType id_node)
{
  checkN2C();
  return m_N2C[m_LNodes[id_node]].size();
}

inline int MeshPartition::n2cGL(vtkIdType id_node, int j)
{
  checkN2C();
  return m_N2C[m_LNodes[id_node]][j];
}

inline vtkIdType MeshPartition::n2cGG(vtkIdType id_node, int j)
{
  checkN2C();
  int i_cell = m_N2C[m_LNodes[id_node]][j];
  if(i_cell<0) return(-1);
  else return m_Cells[i_cell];
}

inline int MeshPartition::c2cLSize(int i_cells)
{
  checkC2C();
  return m_C2C[i_cells].size();
}

inline int MeshPartition::c2cLL(int i_cells, int j)
{
  checkC2C();
  return m_C2C[i_cells][j];
}

inline vtkIdType MeshPartition::c2cLG(int i_cells, int j)
{
  checkC2C();
  int i_cell = m_C2C[i_cells][j];
  if(i_cell<0) return(-1);
  else return m_Cells[i_cell];
}

inline int MeshPartition::c2cGSize(vtkIdType id_cell)
{
  checkC2C();
  checkLCells();
  return m_C2C[m_LCells[id_cell]].size();
}

inline int MeshPartition::c2cGL(vtkIdType id_cell, int j)
{
  checkC2C();
  checkLCells();
  return m_C2C[m_LCells[id_cell]][j];
}

inline vtkIdType MeshPartition::c2cGG(vtkIdType id_cell, int j)
{
  checkC2C();
  checkLCells();
  int i_cell = m_C2C[m_LCells[id_cell]][j];
  if(i_cell<0) return(-1);
  else return m_Cells[i_cell];
}

inline int MeshPartition::getNumberOfCells()
{
  return m_Cells.size();
}

inline int MeshPartition::getNumberOfNodes()
{
  checkNodes();
  return m_Nodes.size();
}

inline int MeshPartition::localNode(vtkIdType id_node)
{
  checkLNodes();
  return m_LNodes[id_node];
}

inline vtkIdType MeshPartition::globalNode(int i)
{
  checkNodes();
  return m_Nodes[i];
}

inline int MeshPartition::localCell(vtkIdType id_cell)
{
  checkLCells();
  return m_LCells[id_cell];
}

inline vtkIdType MeshPartition::globalCell(int i)
{
  if(i<0) return(-1);
  else return m_Cells[i];
}

#endif // MESHPARTITION_H
