//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "guiselectboundarycodes.h"
#include <vtkIntArray.h>
#include <vtkCellData.h>

GuiSelectBoundaryCodes::GuiSelectBoundaryCodes()
{
  disableAutoSet();
};

void GuiSelectBoundaryCodes::setDisplayBoundaryCodes(const QSet<int> &bcs)
{
  display_boundary_codes.clear();
  int bc;
  foreach(bc, bcs) {
    display_boundary_codes.insert(bc);
  };
};

void GuiSelectBoundaryCodes::before()
{
  for (QSet<int>::iterator i = boundary_codes.begin(); i != boundary_codes.end(); ++i) {
    bool checked = display_boundary_codes.contains(*i);
    addListItem(ui.listWidget,*i,checked);
  };
  connect(ui.pushButtonSelect, SIGNAL(clicked()), this, SLOT(selectAll()));
  connect(ui.pushButtonDeselect, SIGNAL(clicked()), this, SLOT(deselectAll()));
};

void GuiSelectBoundaryCodes::operate()
{
  display_boundary_codes.clear();
  for (QSet<int>::iterator i = boundary_codes.begin(); i != boundary_codes.end(); ++i) {
    if (checkListItem(ui.listWidget,*i)) {
      display_boundary_codes.insert(*i);
    };
  };
};

void GuiSelectBoundaryCodes::selectAll()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Checked);
  };
};

void GuiSelectBoundaryCodes::deselectAll()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  };
};

void GuiSelectBoundaryCodes::getSelectedBoundaryCodes(QSet<int> &bcs)
{
  bcs.clear();
  foreach(int bc, display_boundary_codes) {
    bcs.insert(bc);
  };
};


