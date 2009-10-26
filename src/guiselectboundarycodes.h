//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#ifndef guiselectboundarycodes_H
#define guiselectboundarycodes_H

#include "ui_guiselectboundarycodes.h"
#include "dialogoperation.h"

class GuiSelectBoundaryCodes : public DialogOperation<Ui::GuiSelectBoundaryCodes, Operation>
{
  
  Q_OBJECT;
  
private: // attributes
  
  QSet<int> m_DisplayBoundaryCodes;
  
private slots:
  
  void selectAll();
  void deselectAll();
  void saveSelectionAsGrid();
  
protected: // methods
  
  virtual void before();
  virtual void operate();
  
public:
  
  GuiSelectBoundaryCodes();
  void setDisplayBoundaryCodes(const QSet<int> &bcs);
  void getSelectedBoundaryCodes(QSet<int> &bcs);
  
};

#endif
