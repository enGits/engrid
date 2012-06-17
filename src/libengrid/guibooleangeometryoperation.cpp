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

#include "guibooleangeometryoperation.h"
#include "booleangeometryoperation.h"

GuiBooleanGeometryOperation::GuiBooleanGeometryOperation()
{
}

void GuiBooleanGeometryOperation::before()
{
  populateBoundaryCodes(m_Ui.m_ListWidgetBCs1);
  populateBoundaryCodes(m_Ui.m_ListWidgetBCs2);
}

void GuiBooleanGeometryOperation::operate()
{
  QSet<int> bcs1, bcs2;
  getSelectedItems(m_Ui.m_ListWidgetBCs1, bcs1);
  getSelectedItems(m_Ui.m_ListWidgetBCs2, bcs2);
  int s1, s2;
  if      (m_Ui.m_ComboBox->currentText() == "add")      { s1 =  1; s2 =  1; }
  else if (m_Ui.m_ComboBox->currentText() == "subtract") { s1 = -1; s2 =  1; }
  else                                                   { s1 = -1; s2 = -1; }
  BooleanGeometryOperation bool_op(m_Grid, bcs1, bcs2, s1, s2);
  bool_op.setNumCutLayers(m_Ui.m_SpinBoxCutLayers->value());
  bool_op();
}
