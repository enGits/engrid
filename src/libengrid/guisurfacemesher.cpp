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
