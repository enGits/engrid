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

class SurfaceOperation : public Operation
{
public:
  SurfaceOperation();
  virtual void operate();
  
    /**
   * Returns a QVector containing 2 neighbour points to which the point Boss can snap.
   * This is used for removing boundary/feature edge vertices without destroying the geometry.
   */
  bool getNeighbours(vtkIdType Boss, QVector <vtkIdType>& Peons);
  
  ///The same for boundary codes!
  bool getNeighbours_BC(vtkIdType Boss, QVector <vtkIdType>& Peons);
  
  vtkIdType FindSnapPoint(vtkUnstructuredGrid *src, vtkIdType DeadNode,QSet <vtkIdType> & DeadCells,QSet <vtkIdType> & MutatedCells,QSet <vtkIdType> & MutilatedCells, int& N_newpoints, int& N_newcells);
  
  int UpdateCurrentMeshDensity();
  int UpdateNodeType_all();
  int UpdateNodeType();
  
  bool DeletePoint(vtkUnstructuredGrid *src, vtkIdType DeadNode, int& N_newpoints, int& N_newcells);
  bool DeleteSetOfPoints(vtkUnstructuredGrid *src, QSet <vtkIdType> DeadNodes, int& N_newpoints, int& N_newcells);
  int NumberOfCommonPoints(vtkIdType node1, vtkIdType node2, bool& IsTetra);
  
  //--------------------------------------
  //Special for UpdateNodeType_all
  void setConvergence( double C ) { Convergence=C; };
  void setNumberOfIterations( int N ) { NumberOfIterations=N; };
  void setRelaxationFactor( double RF ) { RelaxationFactor=RF; };
  void setFeatureEdgeSmoothing( int FES ) { FeatureEdgeSmoothing=FES; };
  int getFeatureEdgeSmoothing() { return(FeatureEdgeSmoothing); };
  void setFeatureAngle( double FA ) { FeatureAngle=FA; };
  void setEdgeAngle( double EA ) { EdgeAngle=EA; };
  void setBoundarySmoothing( int BS ) { BoundarySmoothing=BS; };
  void setGenerateErrorScalars( int GES ) { GenerateErrorScalars=GES; };
  void setGenerateErrorVectors( int GEV ) { GenerateErrorVectors=GEV; };
  //--------------------------------------
  
  ///Returns the average distance of id_node to its neighbours
  double CurrentVertexAvgDist(vtkIdType id_node);
  
  ///Returns 1/CurrentVertexAvgDist(id_node)
  double CurrentMeshDensity(vtkIdType id_node);
  
  ///Returns the average of 1./node_meshdensity_desired of the neighbours of id_node
  double DesiredVertexAvgDist(vtkIdType id_node);
  
  ///Returns the average of node_meshdensity_desired of the neighbours of id_node
  double DesiredMeshDensity(vtkIdType id_node);
  
  ///Returns the set of boundary codes next to this node
  QSet <int> getBCset(vtkIdType a_node);
  
  ///Returns the node type
  char getNodeType(vtkIdType a_node);
  
  ///Returns the type of the edge [a_node1,a_node2] based on the topology
  char getEdgeType(vtkIdType a_node1, vtkIdType a_node2);
  
  ///Returns the type of the edge [a_node1,a_node2] based on the the type of the two nodes
  char getEdgeType_from_nodes(vtkIdType a_node1, vtkIdType a_node2);
  
  ///passes a vector containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells(vtkIdType id_node1, vtkIdType id_node2,QVector <vtkIdType> &EdgeCells);
  
  ///passes a set containing the cells surrounding edge [id_node1,id_node2] by reference and returns its size
  int getEdgeCells(vtkIdType id_node1, vtkIdType id_node2,QSet <vtkIdType> &EdgeCells);
  
  /// Get VertexMeshDensity object
  VertexMeshDensity getVMD(vtkIdType node);
  
  ///Sets the projection surface and builds the corresponding vtkCellLocator
  void setSource(vtkUnstructuredGrid *a_ProjectionSurface);
  
  ///Sets the vtkCellLocator and ProjectionSurface pointers
  void set_CellLocator_and_ProjectionSurface(vtkCellLocator *a_CellLocator, vtkUnstructuredGrid *a_ProjectionSurface);
  
  /// Projection function
  vec3_t project(vec3_t a_M);
  
  /// Delete CellLocator and ProjectionSurface
  void delete_CellLocator_and_ProjectionSurface();
  
  ///returns the stencil containing id_cell1 and the neighbour cell on side j1 of id_cell1
  stencil_t getStencil(vtkIdType id_cell1, int j1);
  
  ///returns the closest node to a_id_node
  vtkIdType getClosestNode(vtkIdType a_id_node,vtkUnstructuredGrid* a_grid);
  
  ///returns the farthest node from a_id_node
  vtkIdType getFarthestNode(vtkIdType a_id_node,vtkUnstructuredGrid* a_grid);
  
  //---------------------------------------------------
  //Utility functions used in Roland's formulas
  //Should be renamed to be more explicit
  //Some could be moved into geometrytools
  //Some are pretty useless
  
  ///perimeter
  double Um(vtkIdType D);
  /// area of the circumscribed circle of the triangle
  double A_U(vtkIdType D);
  /// triangle area
  double A_D(vtkIdType D);
  /// triangle neighbours
  double DN(int i,vtkIdType D);
  /// number of edges
  double nk(vtkIdType P);
  
  double G_k(vtkIdType node);
  /// triangle nodes
  double DK(int i,vtkIdType D);
  
  vtkIdType KK(int i,vtkIdType j,vtkIdType K);
  
  double L_k(vtkIdType j,vtkIdType K);// node1 K, node2 j
  
  double Q_L(vtkIdType D);
  
  double Q_L1(vtkIdType P);
  
  double Q_L2(vtkIdType P);
  
  double T_min(int w);
//---------------------------------------------------
};

#endif
