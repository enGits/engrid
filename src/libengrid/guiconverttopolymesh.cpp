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
#include "guiconverttopolymesh.h"
#include "converttopolymesh.h"

void GuiConvertToPolyMesh::operate()
{
  ConvertToPolyMesh convert;
  if (m_Ui.checkBoxOptimise->isChecked()) {
    convert.setOptimiseOn();
  }
  if (m_Ui.checkBoxSplitConcaveFaces->isChecked()) {
    convert.setSplitFacesOn();
  }
  if (m_Ui.checkBoxSplitConcaveCells->isChecked()) {
    convert.setSplitCellsOn();
  }
  convert.setPullInFactor(m_Ui.horizontalSliderPullIn->value()*0.01);
  convert();
}
