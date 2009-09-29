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
#ifndef FIXCADGEOMETRY_H
#define FIXCADGEOMETRY_H

#include "surfacealgorithm.h"

class fixCadGeometry: public SurfaceAlgorithm
{
  
protected: // methods
  
  virtual void operate();
  
  
public: // methods
  
  fixCadGeometry();
  void mesher();
  void setDesiredLength(double L=9000);
  void customUpdateNodeInfo(bool update_type = false);
  
  /// Returns the node type
  char custom_getNodeType( vtkIdType a_node, bool allow_feature_edge_vertices = false, bool fix_unselected = true );
  
  /// Returns the type of the edge [a_node1,a_node2] based on the topology
  char custom_getEdgeType( vtkIdType a_node1, vtkIdType a_node2, bool allow_feature_edge_vertices, bool fix_unselected );
  
};

#endif
