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

#include "guicreatehexcore.h"

GuiCreateHexCore::GuiCreateHexCore()
{
}

void GuiCreateHexCore::before()
{
  vec3_t x1( 1e99,  1e99,  1e99);
  vec3_t x2(-1e99, -1e99, -1e99);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    for (int i = 0; i < 3; ++i) {
      x1[i] = min(x1[i], x[i]);
      x2[i] = max(x2[i], x[i]);
    }
  }
  double xmin = min(x1[0], min(x1[1], x1[2]));
  double xmax = max(x2[0], max(x2[1], x2[2]));
  m_X1 = vec3_t(xmin, xmin, xmin);
  m_X2 = vec3_t(xmax, xmax, xmax);
  QString num;
  vec3_t xi = 0.5*(x1 + x2);
  num.setNum(xi[0]);   m_Ui.lineEditCiX->setText(num);
  num.setNum(xi[1]);   m_Ui.lineEditCiY->setText(num);
  num.setNum(xi[2]);   m_Ui.lineEditCiZ->setText(num);
  num.setNum(m_X1[0]); m_Ui.lineEditX1->setText(num);
  num.setNum(m_X1[1]); m_Ui.lineEditY1->setText(num);
  num.setNum(m_X1[2]); m_Ui.lineEditZ1->setText(num);
  num.setNum(m_X2[0]); m_Ui.lineEditX2->setText(num);
  num.setNum(m_X2[1]); m_Ui.lineEditY2->setText(num);
  num.setNum(m_X2[2]); m_Ui.lineEditZ2->setText(num);
}

void GuiCreateHexCore::operate()
{
  vec3_t xi(m_Ui.lineEditCiX->text().toDouble(), m_Ui.lineEditCiY->text().toDouble(), m_Ui.lineEditCiZ->text().toDouble());
  m_X1 = vec3_t(m_Ui.lineEditX1->text().toDouble(), m_Ui.lineEditY1->text().toDouble(), m_Ui.lineEditZ1->text().toDouble());
  m_X2 = vec3_t(m_Ui.lineEditX2->text().toDouble(), m_Ui.lineEditY2->text().toDouble(), m_Ui.lineEditZ2->text().toDouble());
  CreateHexCore create_hex_core(m_X1, m_X2, xi);
  create_hex_core();
}

