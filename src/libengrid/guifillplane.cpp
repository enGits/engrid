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

#include "guifillplane.h"
#include "fillplane.h"

void GuiFillPlane::operate()
{
  vec3_t x, n;
  x[0] = m_Ui.m_LineEditX->text().toDouble();
  x[1] = m_Ui.m_LineEditY->text().toDouble();
  x[2] = m_Ui.m_LineEditZ->text().toDouble();
  n[0] = m_Ui.m_LineEditNx->text().toDouble();
  n[1] = m_Ui.m_LineEditNy->text().toDouble();
  n[2] = m_Ui.m_LineEditNz->text().toDouble();
  double tol = m_Ui.m_LineEditPrec->text().toDouble();
  n.normalise();
  FillPlane fill;
  fill.setOrigin(x);
  fill.setNormal(n);
  fill.setTolerance(tol);
  fill();
}
