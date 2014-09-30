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
// 
#ifndef GEOMETRYSMOOTHER_H
#define GEOMETRYSMOOTHER_H

#include "operation.h"

class GeometrySmoother;


class GeometrySmoother : public Operation
{

protected: // attributes

  double m_RelaxationFactor;
  double m_CornerAngle;
  int    m_NumIterations;


protected: //  methods

  virtual void operate();


public: // methods

  GeometrySmoother();

  void setRelaxationFactor(double relaxation_factor) { m_RelaxationFactor = relaxation_factor; }
  void setCornerAngle(double edge_angle)             { m_CornerAngle = edge_angle; }
  void setNumberOfIterations(int num)                { m_NumIterations = num; }


};

#endif // GEOMETRYSMOOTHER_H
