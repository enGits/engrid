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
#ifndef SURFACEMESHSMOOTHER_H
#define SURFACEMESHSMOOTHER_H

class SurfaceMeshSmoother;

#include "optimisation.h"
#include "operation.h"
#include "cadinterface.h"

class SurfaceMeshSmoother : public Operation
{

private: // attributes

  //vtkIdType                    m_IdNode;
  bool                         m_UseEstimatedPlane;
  vec3_t                       m_X0;
  vec3_t                       m_N0;
  CadInterface                *m_Cad;
  QList<QPair<vec2_t,vec2_t> > m_Limit;
  mat3_t                       m_M32;
  mat3_t                       m_M23;
  double                       m_MinLength;
  double                       m_Prec;
  int                          m_NumSteps;
  bool                         m_UseSimpleCentreScheme;


protected: // methods

  vec2_t transform32(vec3_t x);
  vec3_t transform23(vec2_t x);

  virtual void operate();
  virtual double error(vec2_t x);

  vec2_t gradError(vec2_t x);
  void computeLimits(vec2_t x, vec2_t v, double &k1, double &k2);
  vec2_t iteration(vec2_t x);


public:

  SurfaceMeshSmoother();

  void prepareEstimatedPlane();
  void prepareCadInterface(CadInterface *cad);
  vec3_t smoothNode(vtkIdType id_node);
  void useSimpleCentreScheme() { m_UseSimpleCentreScheme = true; }

};

#endif // SURFACEMESHSMOOTHER_H
