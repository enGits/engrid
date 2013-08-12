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
#ifndef CONVERTTOPOLYMESH_H
#define CONVERTTOPOLYMESH_H

#include "operation.h"

class ConvertToPolyMesh : public Operation
{

protected: // attributes

  bool   m_Optimise;
  bool   m_SplitFaces;
  bool   m_SplitCells;
  double m_PullInFactor;


protected: // methods

  virtual void operate();


public: // methods;

  ConvertToPolyMesh();

  void setOptimiseOn()    { m_Optimise   = true; }
  void setOptimiseOff()   { m_Optimise   = false; }
  void setSplitFacesOn()  { m_SplitFaces = true; }
  void setSplitFacesOff() { m_SplitFaces = false; }
  void setSplitCellsOn()  { m_SplitCells = true; }
  void setSplitCellsOff() { m_SplitCells = false; }

  void setPullInFactor(double f) { m_PullInFactor = f; }


};

#endif // CONVERTTOPOLYMESH_H
