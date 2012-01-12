// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#ifndef setboundarycode_H
#define setboundarycode_H

class SetBoundaryCode;

#include "cellneighbouriterator.h"

class SetBoundaryCode : public CellNeighbourIterator
{
  
private: // attributes
  
  double feature_angle;
  int    boundary_code;
  bool   ProcessAll;
  bool   SelectAllVisible;
  bool   OnlyPickedCell;
  bool   OnlyPickedCellAndNeighbours;
  
protected: // methods
  
  virtual void pass1();
  virtual void pass2();
  
public: // methods
  
  SetBoundaryCode();
  void setFeatureAngle(double fa) { feature_angle = fa; }
  void setBC(int bc) { boundary_code = bc; }
  void setProcessAll(bool b) { ProcessAll=b; }
  void setSelectAllVisible(bool b) { SelectAllVisible=b; }
  void setOnlyPickedCell(bool b) { OnlyPickedCell=b; }
  void setOnlyPickedCellAndNeighbours(bool b) { OnlyPickedCellAndNeighbours=b; }
  void set(bool b) { OnlyPickedCell=b; }
  
};

#endif
