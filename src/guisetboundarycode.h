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
#ifndef guisetboundarycode_H
#define guisetboundarycode_H

#include <QButtonGroup>
#include <QRadioButton>

class GuisetBoundaryCode;

#include "dialogoperation.h"
#include "ui_guisetboundarycode.h"

class GuisetBoundaryCode : public DialogOperation<Ui::GuisetBoundaryCode, Operation>
{
  
  Q_OBJECT;
private:
  QButtonGroup* buttongroup;
  QRadioButton* radioButton_ProcessOnlyVisible;
  QRadioButton* radioButton_ProcessAll;
  QRadioButton* radioButton_SelectAllVisible;
  QRadioButton* radioButton_OnlyPickedCell;
  QRadioButton* radioButton_OnlyPickedCellAndNeighbours;
  
protected: // methods
  
  virtual void before();
  virtual void operate();
  
};

#endif
