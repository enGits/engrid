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

#ifndef STITCHHOLES_H
#define STITCHHOLES_H

#include "fillplane.h"
#include "cadinterface.h"

class StitchHoles : public FillPlane
{

private: // attributes

  int           m_Bc;
  CadInterface *m_Cad;

  QList<vec3_t> m_X2;
  QList<vec3_t> m_X3;

  int m_Counter;

protected: // methods

  QList<vtkIdType> getNextHole();
  void stitchHole(QList<vtkIdType> loop_nodes);
  vec3_t transformFromPlane(vec3_t x);

  virtual void operate();


public:

  StitchHoles(int bc);

};

#endif // STITCHHOLES_H
