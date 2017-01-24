// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2016 enGits GmbH                                      +
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

#ifndef PADSURFACE_H
#define PADSURFACE_H

#include "operation.h"

class PadSurface : public Operation
{

protected: // attributes

  QSet<int> m_BCs;
  int       m_NewBC;
  bool      m_Relative;
  double    m_Distance;


protected: // methods

  virtual void operate();


public:

  PadSurface();
  void addBC(int bc) { m_BCs.insert(bc); }
  void addBC(BoundaryCondition bc) { m_BCs.insert(bc.getCode()); }
  void setNewBC(int bc) { m_NewBC = bc; }
  void setDistance(double d) { m_Distance = d; }
  void relativeOn() { m_Relative = true; }
  void relativeOff() { m_Relative = false; }

};

#endif // PADSURFACE_H
