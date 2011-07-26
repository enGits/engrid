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
#include "guisurfacemesher.h"

void GuiSurfaceMesher::before()
{
  ui.m_SpinBoxIterations->setValue(m_NumMaxIter);
  ui.m_SpinBoxDelaunaySweeps->setValue(m_NumDelaunaySweeps);
  ui.m_SpinBoxSmoothSteps->setValue(m_NumSmoothSteps);
  ui.m_CheckBoxSmooth->setChecked(m_CorrectCurvature);
}

void GuiSurfaceMesher::operate()
{
  setMaxNumIterations(ui.m_SpinBoxIterations->value());
  setNumDelaunaySweeps(ui.m_SpinBoxDelaunaySweeps->value());
  setNumSmoothSteps(ui.m_SpinBoxSmoothSteps->value());
  setCorrectCurvature(ui.m_CheckBoxSmooth->isChecked());
  SurfaceMesher::operate();
}
