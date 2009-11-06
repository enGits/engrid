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
#ifndef SURFACEOPERATION_H
#define SURFACEOPERATION_H

#include <operation.h>

//==============================================
/// Special structure for working on two linked cells
// struct stencil_t {
//   vtkIdType id_cell1;
//   vtkIdType id_cell2;
//   vtkIdType p[4];
//   bool sameBC;         ///< do both cells have the same BCs?
//   bool twocells;       ///< Do we have 2 cells?
//   char neighbour_type; ///< What's the type of the neighbour cell?
// };

//==============================================
/** Special structure for working on two (or more) linked cells
 * WARNING: Only one exterior node is stored for non-triangular cells!
 */
struct stencil_t
{
  QVector<vtkIdType> id_cell;   ///< stencil cells
  QVector<vtkIdType> id_node;   ///< "exterior nodes" of the stencil cells
  QVector<vtkIdType> type_cell; ///< type of the stencil cells
  vtkIdType p1, p2;             ///< points forming the central edge of the stencil
  bool sameBC;                  ///< do all cells have the same BCs?
};

/// Prints out stencil information
ostream& operator<<( ostream &out, stencil_t S );

//==============================================

class SurfaceOperation : public Operation
{

private:

  ///Vector used to store the "Potential Snap Points" of each point, i.e. the neighbour points belonging to the same edge (boundary or feature) in case of edge points, all neighbour points in case of simple points and the points belonging to edges in case of fixed points
  QVector < QVector <vtkIdType> > m_PotentialSnapPoints;


protected:

  ///\todo Remove useless attributes
  //attributes for determining node types and for smoothing operations
  double m_Convergence;
  int    m_NumberOfIterations;
  double m_RelaxationFactor;
  //  int    m_AllowFeatureEdgeVertices; ///< if set to 0, feature edge vertices will be deactivated. Use setm_AllowFeatureEdgeVertices(int) to set it.
  double m_FeatureAngle;
  double m_EdgeAngle;
  int    m_BoundarySmoothing;


public:

  SurfaceOperation();
  virtual void operate();

  QVector <vtkIdType> getPotentialSnapPoints( vtkIdType id_node ); ///< Returns a QVector containing neighbour points to which the point id_node can snap.

  int UpdateCurrentMeshDensity();

  /// Updates the m_PotentialSnapPoints structure + updates node types if desired (faster than loop through nodes with getNodeType)
  int UpdatePotentialSnapPoints(bool update_node_types, bool fix_unselected = true);

  //--------------------------------------
  //Special for UpdatePotentialSnapPoints
  void setConvergence( double C )           { m_Convergence = C; }
  void setNumberOfIterations( int N )       { m_NumberOfIterations = N; }
  void setRelaxationFactor( double RF )     { m_RelaxationFactor = RF; }
  //void setAllowFeatureEdgeVertices( int x ) {  = x; } ///< If x = 0, feature edge vertices will be deactivated.
  //int getAllowFeatureEdgeVertices() { return( m_AllowFeatureEdgeVertices ); }
  void setFeatureAngle(double FA)   { m_FeatureAngle = FA; }
  void setEdgeAngle(double EA)      { m_EdgeAngle = EA; }
  void setBoundarySmoothing(int BS) { m_BoundarySmoothing = BS; }
  //--------------------------------------

  /// Returns the average distance of id_node to its neighbours
  double CurrentVertexAvgDist( vtkIdType id_node );

  /// Returns 1/CurrentVertexAvgDist(id_node)
  double CurrentMeshDensity( vtkIdType id_node );

  /// Returns the average of 1./node_meshdensity_desired of the neighbours of id_node
  double DesiredVertexAvgDist( vtkIdType id_node );

  /// Returns the average of node_meshdensity_desired of the neighbours of id_node
  double DesiredMeshDensity( vtkIdType id_node );

  /// Returns the set of boundary codes next to this node
  QSet <int> getBCset( vtkIdType a_node );

  /// Returns the node type
  char getNodeType(vtkIdType a_node, bool fix_unselected = true);

  /// Returns the type of the edge [a_node1,a_node2] based on the topology
  char getEdgeType(vtkIdType a_node1, vtkIdType a_node2, bool fix_unselected);

  // Returns the type of the edge [a_node1,a_node2] based on the the type of the two nodes
  // deprecated?
//   char getEdgeType_from_nodes( vtkIdType a_node1, vtkIdType a_node2 );

  /// passes a vector containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QVector <vtkIdType> &EdgeCells );

  /// passes a set containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QSet <vtkIdType> &EdgeCells );

  /// Get VertexMeshDensity object
  VertexMeshDensity getVMD( vtkIdType node );

  /// returns the stencil containing id_cell1 and the neighbour cell on side j1 of id_cell1
  stencil_t getStencil( vtkIdType id_cell1, int j1 );

  /// returns the closest neighbour node of id_node
  vtkIdType getClosestNode( vtkIdType id_node );

  /// returns the farthest neighbour node of id_node
  vtkIdType getFarthestNode( vtkIdType id_node );

  //---------------------------------------------------
  //Utility functions used in Roland's formulas

  /// perimeter
  double perimeter( vtkIdType id_cell );

  /// desired edge length for id_node
  double desiredEdgeLength( vtkIdType id_node );

  /// mean desired edge length for id_cell
  double meanDesiredEdgeLength( vtkIdType id_cell );

  /// perimeter / sum of the desired edge lengths
  double Q_L( vtkIdType id_cell );

  /// sum(2*edgelength,edges(id_node))/sum(desired edgelengths of each edgepoint,edges(id_node))
  double Q_L1( vtkIdType id_node );

  /// minimum of sum(2*edgelength)/sum(desired edgelengths of each edgepoint) for each edge of id_node
  double Q_L2( vtkIdType id_node );

  /// Value to minimize for mesh smoothing. w allows putting more weight on the form or the area of triangles.
  double T_min( int w );

  //---------------------------------------------------
};

#endif
