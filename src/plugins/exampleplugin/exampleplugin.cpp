//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "exampleplugin.h"

#include <QtGui>

void ExamplePlugin::before()
{
  ui.lineEdit_Scaling_X->setText("1");
  ui.lineEdit_Scaling_Y->setText("1");
  ui.lineEdit_Scaling_Z->setText("1");
}

void ExamplePlugin::operate()
{
  double sx = ui.lineEdit_Scaling_X->text().toDouble();
  double sy = ui.lineEdit_Scaling_Y->text().toDouble();
  double sz = ui.lineEdit_Scaling_Z->text().toDouble();
  l2g_t nodes = getPartNodes();
  foreach(vtkIdType id_node, nodes) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    x[0] *= sx;
    x[1] *= sy;
    x[2] *= sz;    
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
  m_Grid->Modified();
};

Q_EXPORT_PLUGIN2(exampleplugin, ExamplePlugin)
