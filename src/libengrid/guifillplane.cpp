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

#include "guifillplane.h"
#include "fillplane.h"
#include "boundarycondition.h"
#include "guimainwindow.h"

void GuiFillPlane::fillPlane(vec3_t x, vec3_t n, bool inverse, QString name)
{
  double tol = m_Ui.m_LineEditDistTol->text().toDouble();
  n.normalise();
  FillPlane fill;
  fill.setOrigin(x);
  fill.setNormal(n);
  fill.setDistanceTolerance(m_Ui.m_LineEditDistTol->text().toDouble());
  fill.setAngularTolerance(GeometryTools::deg2rad(m_Ui.m_LineEditAngTol->text().toDouble()));
  if (inverse) {
    fill.setInverseDirectionOn();
  } else {
    fill.setInverseDirectionOff();
  }
  fill();
  BoundaryCondition bc(name, "patch", fill.getBC());
  GuiMainWindow::pointer()->setBC(fill.getBC(), bc);
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
}

void GuiFillPlane::operate()
{
  if (m_Ui.m_CheckBoxCustom->isChecked()) {
    vec3_t x, n;
    x[0] = m_Ui.m_LineEditX->text().toDouble();
    x[1] = m_Ui.m_LineEditY->text().toDouble();
    x[2] = m_Ui.m_LineEditZ->text().toDouble();
    n[0] = m_Ui.m_LineEditNx->text().toDouble();
    n[1] = m_Ui.m_LineEditNy->text().toDouble();
    n[2] = m_Ui.m_LineEditNz->text().toDouble();
    fillPlane(x, n, false, "filled_plane");
  }
  if (m_Ui.m_CheckBoxXY->isChecked()) {
    vec3_t x(0, 0, 0);
    vec3_t n(0, 0, 1);
    fillPlane(x, n, m_Ui.m_CheckBoxXYInv->isChecked(), "XY_plane");
  }
  if (m_Ui.m_CheckBoxYZ->isChecked()) {
    vec3_t x(0, 0, 0);
    vec3_t n(1, 0, 0);
    fillPlane(x, n, m_Ui.m_CheckBoxYZInv->isChecked(), "YZ_plane");
  }
  if (m_Ui.m_CheckBoxZX->isChecked()) {
    vec3_t x(0, 0, 0);
    vec3_t n(0, 1, 0);
    fillPlane(x, n, m_Ui.m_CheckBoxZXInv->isChecked(), "ZX_plane");
  }
}
