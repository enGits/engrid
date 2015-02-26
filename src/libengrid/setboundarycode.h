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
#ifndef setboundarycode_H
#define setboundarycode_H

class SetBoundaryCode;

#include "cellneighbouriterator.h"

class SetBoundaryCode : public CellNeighbourIterator
{
  
private: // attributes
  
  double m_FeatureAngle;
  int    m_NewBoundaryCode;
  int    m_OldBoundaryCode;
  bool   m_ProcessAll;
  bool   m_SelectAllVisible;
  bool   m_OnlyPickedCell;
  bool   m_OnlyPickedCellAndNeighbours;
  
protected: // methods
  
  virtual void pass1();
  virtual void pass2();
  
public: // methods
  
  SetBoundaryCode();
  void setFeatureAngle(double fa) { m_FeatureAngle = fa; }
  void setNewBC(int bc) { m_NewBoundaryCode = bc; }
  void setOLdBC(int bc) { m_OldBoundaryCode = bc; }
  void setProcessAll(bool b) { m_ProcessAll=b; }
  void setSelectAllVisible(bool b) { m_SelectAllVisible=b; }
  void setOnlyPickedCell(bool b) { m_OnlyPickedCell=b; }
  void setOnlyPickedCellAndNeighbours(bool b) { m_OnlyPickedCellAndNeighbours=b; }

};

#endif
