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
#ifndef SURFACEOPERATION_H
#define SURFACEOPERATION_H

#include <operation.h>

//==============================================

/// Special structure for working on two (or more) linked cells
struct stencil_t
{
  QVector<vtkIdType> id_cell;
  QVector<vtkIdType> id_node;
  QVector<vtkIdType> type_cell;
  vtkIdType p1, p2;
  bool sameBC;         ///< do all cells have the same BCs?
};

/// Prints out stencil information
ostream& operator<<( ostream &out, stencil_t S );

//==============================================

class SurfaceOperation : public Operation
{

private:

  ///Vector used to store the "Potential Snap Points" of each point, i.e. the neighbour points belonging to the same edge (boundary or feature) in case of edge points, all neighbour points in case of simple points and the points belonging to edges in case of fixed points
  QVector < QVector <vtkIdType> > m_PotentialSnapPoints;


protected: // attributes

  double m_FeatureAngle;
  double m_EdgeAngle;
  int    m_BoundarySmoothing;

  QVector<vec3_t> m_NodeNormal; ///< node normal vectors
  double m_StretchingFactor;


protected: // methods

  void computeNormals();
  double normalIrregularity(vtkIdType id_node);

public:

  SurfaceOperation();
  virtual void operate();

  QVector <vtkIdType> getPotentialSnapPoints( vtkIdType id_node ); ///< Returns a QVector containing neighbour points to which the point id_node can snap.

  int UpdateCurrentMeshDensity();

  /// Updates the m_PotentialSnapPoints structure + updates node types if desired (faster than loop through nodes with getNodeType)
  int UpdatePotentialSnapPoints(bool update_node_types, bool fix_unselected = true);

  void setFeatureAngle(double FA)   { m_FeatureAngle = FA; }
  void setEdgeAngle(double EA)      { m_EdgeAngle = EA; }
  void setBoundarySmoothing(int BS) { m_BoundarySmoothing = BS; }


  double currentVertexAvgDist(vtkIdType id_node);                 ///< Returns the average distance of id_node to its neighbours
  double CurrentMeshDensity( vtkIdType id_node );                 ///< Returns 1/CurrentVertexAvgDist(id_node)
  char getNodeType(vtkIdType a_node, bool fix_unselected = true); ///< Returns the node type

  /**
   * Get the type of a given edge based on the topology.
   * @param id_node1 first node of the edge
   * @param id_node2 second node of the edge
   * @param fix_unselected fix all edges which belong to unselected boundary codes
   * @return the type of the edge
   */
  char getEdgeType(vtkIdType id_node1, vtkIdType od_node2, bool fix_unselected);

  /// passes a vector containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QVector <vtkIdType> &EdgeCells );

  /// passes a set containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QSet <vtkIdType> &EdgeCells );

  /// Get VertexMeshDensity object
  VertexMeshDensity getVMD(vtkIdType id_node);

  /// returns the stencil containing id_cell1 and the neighbour cell on side j1 of id_cell1
  stencil_t getStencil(vtkIdType id_cell1, int j1);

  /// desired edge length for id_node
  double desiredEdgeLength(vtkIdType id_node);

  /// mean desired edge length for id_cell
  double meanDesiredEdgeLength(vtkIdType id_cell);

  bool isCell(vtkIdType id_node1, vtkIdType id_node2, vtkIdType id_node3);

  void setStretchingFactor(double sf) { m_StretchingFactor = sf; }

};

#endif
