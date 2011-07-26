// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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

class GuiSetBoundaryCode;

#include "dialogoperation.h"
#include "ui_guisetboundarycode.h"

class GuiSetBoundaryCode : public DialogOperation<Ui::GuiSetBoundaryCode, Operation>
{
  
  Q_OBJECT;

private:

  QButtonGroup* m_ButtonGroup;
  QRadioButton* m_RadioButtonProcessAll;
  QRadioButton* m_RadioButtonProcessOnlyVisible;
  QRadioButton* m_RadioButtonSelectAllVisible;
  QRadioButton* m_RadioButtonOnlyPickedCell;
  QRadioButton* m_RadioButtonOnlyPickedCellAndNeighbours;
  QRadioButton* m_RadioButtonAuto;

protected: // methods
  
  virtual void before();
  virtual void operate();
  
};

#endif
