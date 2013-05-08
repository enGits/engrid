//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                     +
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
#ifndef SURFACEPROJECTION_H
#define SURFACEPROJECTION_H

#include "surfacealgorithm.h"
#include "cadinterface.h"

class SurfaceProjection : public SurfaceAlgorithm
{

private: // attributes

  vtkUnstructuredGrid* m_FGrid;        ///< the foreground grid to project
  MeshPartition        m_FPart;        ///< MeshPartition for the foreground grid
  CadInterface*        m_CadInterface; ///< the CAD interface providing the geometry description

  vec3_t m_LastNormal;
  double m_LastRadius;
  bool   m_ForceRay;
  bool   m_Failed;

public: // methods

  SurfaceProjection(CadInterface* cad_interface);

  void    setForegroundGrid(vtkUnstructuredGrid* grid);

  virtual vec3_t projectNode(vec3_t x, vtkIdType id_node = -1, bool correct_curvature = false, vec3_t v = vec3_t(0,0,0),
                             bool strict_direction = false, bool allow_search = true);
  virtual vec3_t snapNode(vec3_t x, vtkIdType id_node = -1, bool correct_curvature = false);

  double    getRadius(vtkIdType id_node);
  vec3_t    lastProjNormal() { return m_LastNormal; }
  double    lastProjRadius() { return m_LastRadius; }
  vec3_t    correctCurvature(vtkIdType, vec3_t x) { return x; }
  vtkIdType lastProjTriangle() { return -1; }
  bool      lastProjFailed() { return m_Failed; }

};

#endif // SURFACEPROJECTION_H
