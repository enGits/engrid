// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#ifndef CREATEBOUNDARYLAYER_H
#define CREATEBOUNDARYLAYER_H

#include "boundarylayeroperation.h"
#include "guimainwindow.h"

class CreateBoundaryLayerShell : public BoundaryLayerOperation
{
private: // attributes

  QVector<vtkIdType> layer_cells;

  int     m_NumIterations;
  bool    m_RemovePoints;
  double  m_RelativeHeight;
  double  m_AbsoluteHeight;
  double  m_Blending;
  QString m_VolumeName;

  VolumeDefinition                     m_VolDef;
  vtkSmartPointer<vtkUnstructuredGrid> m_RestGrid;
  vtkSmartPointer<vtkUnstructuredGrid> m_OriginalGrid;
  vtkSmartPointer<vtkUnstructuredGrid> m_PrismaticGrid;


  /// Boundary codes of the surface we want to remove points on. Normally the one next to the prismatic boundary layer.
  QSet<int> m_LayerAdjacentBoundaryCodes;

  QVector<vtkIdType> m_ShellNodeMap;
  MeshPartition      m_ShellPart;

  QMap<int, vec3_t> m_LayerAdjacentOrigins;
  QMap<int, vec3_t> m_LayerAdjacentNormals;


private: // methods

  void correctAdjacentBC(int bc, vtkUnstructuredGrid *grid);
  void prepare();
  void createLayerNodes(vtkIdType id_node);
  void createPrismaticGrid();


protected: // methods

  virtual void operate();

  void reduceSurface();
  void smoothSurface();


public: // methods

  CreateBoundaryLayerShell();

  void setVolumeName(QString name) { m_VolumeName = name; }
  vtkUnstructuredGrid* getPrismaticGrid() { return m_PrismaticGrid; }
  QVector<int> getBoundaryLayerCodes() { return m_BoundaryLayerCodes; }

};

#endif // CREATEBOUNDARYLAYER_H
