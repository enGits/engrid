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

#include "guicreatehexcore.h"

GuiCreateHexCore::GuiCreateHexCore()
{
}

void GuiCreateHexCore::before()
{

}

void GuiCreateHexCore::operate()
{
  vec3_t x1(ui.lineEditC1X->text().toDouble(), ui.lineEditC1Z->text().toDouble(), ui.lineEditC1Y->text().toDouble());
  vec3_t x2(ui.lineEditC2X->text().toDouble(), ui.lineEditC2Z->text().toDouble(), ui.lineEditC2Y->text().toDouble());
  vec3_t xi(ui.lineEditCiX->text().toDouble(), ui.lineEditCiZ->text().toDouble(), ui.lineEditCiY->text().toDouble());
  CreateHexCore create_hex_core(x1, x2, xi);
  create_hex_core();
}

