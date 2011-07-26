// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#ifndef REMOVEPOINTS_H
#define REMOVEPOINTS_H

#include "surfaceoperation.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCharArray.h>

#include <QVector>
#include <QString>
#include <QTextStream>
#include <QTime>

class RemovePoints : public SurfaceOperation
{

protected:

  int    m_NumRemoved;
  double m_Threshold;
  bool   m_ProtectFeatureEdges;
  bool   m_PerformGeometricChecks;
  bool   m_UpdatePSP;

  QVector<bool> m_IsFeatureNode;
  QVector<bool> m_Fixed;

public:

  RemovePoints();

  virtual void operate();

  int getNumRemoved() { return m_NumRemoved; }
  void setProtectFeatureEdgesOn()  { m_ProtectFeatureEdges = true; }
  void setProtectFeatureEdgesOff() { m_ProtectFeatureEdges = false; }
  void setPerformGeometricChecksOn()  { m_PerformGeometricChecks = true; }
  void setPerformGeometricChecksOff() { m_PerformGeometricChecks = false; }
  void setUpdatePSPOn()  { m_UpdatePSP = true; }
  void setUpdatePSPOff() { m_UpdatePSP = false; }
  void setThreshold(double v) { m_Threshold = v; }
  void fixNodes(const QVector<bool> &fixnodes);

protected:

  void markFeatureEdges();

  /// deletes set of points DeadNodes
  bool DeleteSetOfPoints(const QVector<vtkIdType>& deadnode_vector,
                         const QVector<vtkIdType>& snappoint_vector,
                         const QVector<vtkIdType>& all_deadcells,
                         const QVector<vtkIdType>& all_mutatedcells);
  
  /// returns a valid potential snappoint (checks for flipped cells, etc). If none is found, returns -1.
  vtkIdType FindSnapPoint( vtkIdType DeadNode,
                           QVector<vtkIdType>& DeadCells,
                           QVector<vtkIdType>& MutatedCells,
                           int& N_newpoints, int & N_newcells,
                           const QVector<bool>& marked_nodes);
  
  /// returns true if moving id_node to position P leads to flipped cells
  bool flippedCell(vtkIdType id_node, vec3_t x_new, vtkIdType id_cell);
  bool flippedCell2(vtkIdType id_node, vec3_t x_new);
  
  /// returns number of common neighbour nodes of id_node1 and id_node2. IsTetra becomes true if id_node1 and id_node2 belong to the edge of a tetrahedron.
  int NumberOfCommonPoints( vtkIdType id_node1, vtkIdType id_node2, bool& IsTetra );

  /// returns number of common neighbour nodes of id_node1 and id_node2. IsTetra becomes true if id_node1 and id_node2 belong to the edge of a tetrahedron.
  bool checkForDestroyedVolumes( vtkIdType id_node1, vtkIdType id_node2, int& N_common_points );
  
};

#endif
