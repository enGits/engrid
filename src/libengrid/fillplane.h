//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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

#ifndef FILLPLANE_H
#define FILLPLANE_H

#include "surfaceoperation.h"

class FillPlane : public SurfaceOperation
{

protected: // attributes

  vec3_t m_X0;
  vec3_t m_N;
  vec3_t m_G1;
  vec3_t m_G2;
  double m_DistTol;
  double m_AngleTol;
  bool   m_InverseDirection;
  int    m_BC;

  QVector<vtkIdType> m_NodeMap;


protected: // methods

  void   createEdgesOnPlane(vtkUnstructuredGrid *edge_grid);
  void   closeLoops(vtkUnstructuredGrid *edge_grid);
  vec3_t toPlane(vec3_t x);
  vec3_t fromPlane(vec3_t x);
  bool   isWithinTolerance(vec3_t x);
  void   gridToPlane(vtkUnstructuredGrid *edge_grid);
  void   gridFromPlane(vtkUnstructuredGrid *edge_grid);
  void   triangulate(vtkPolyData *edge_pdata, vtkUnstructuredGrid *tri_grid);
  void   order(vtkUnstructuredGrid *edge_grid, vtkPolyData *edge_pdata);
  void   orderGeometrically(vtkUnstructuredGrid *edge_grid, QList<vtkIdType> &poly_nodes);

  virtual void operate();

public:

  FillPlane();

  void setOrigin(vec3_t x) { m_X0 = x; }
  void setNormal(vec3_t n) { m_N = n; }
  void setDistanceTolerance(double t) { m_DistTol = t; }
  void setAngularTolerance(double a) { m_AngleTol = a; }
  void setInverseDirectionOn() { m_InverseDirection = true; }
  void setInverseDirectionOff() { m_InverseDirection = false; }
  int  getBC() { return m_BC; }

};

#endif // FILLPLANE_H
