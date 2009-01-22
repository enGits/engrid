//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "guieditboundaryconditions.h"

GuiEditBoundaryConditions::GuiEditBoundaryConditions()
{
  bcmap = NULL;
};

void GuiEditBoundaryConditions::before()
{
  if (!bcmap) EG_BUG;
  while (ui.T->rowCount()) ui.T->removeRow(0);
  foreach (int i, boundary_codes) {
    BoundaryCondition bc = (*bcmap)[i];
    ui.T->insertRow(ui.T->rowCount());
    int r = ui.T->rowCount()-1;
    ui.T->setItem(r,0,new QTableWidgetItem());
    ui.T->item(r,0)->setFlags(ui.T->item(r,0)->flags() & (~Qt::ItemIsSelectable));
    ui.T->item(r,0)->setFlags(ui.T->item(r,0)->flags() & (~Qt::ItemIsEditable));
    ui.T->setItem(r,1,new QTableWidgetItem());
    ui.T->setItem(r,2,new QTableWidgetItem());
    QString idx;
    idx.setNum(i);
    ui.T->item(r,0)->setText(idx);
    QString name = bc.getName();
    if (name == "unknown") name = QString("BC") + idx;
    ui.T->item(r,1)->setText(name);
    ui.T->item(r,2)->setText(bc.getType());
  };
  
};

void GuiEditBoundaryConditions::operate()
{
  for (int i = 0; i < ui.T->rowCount(); ++i) {
    BoundaryCondition bc(ui.T->item(i,1)->text(),ui.T->item(i,2)->text());
    (*bcmap)[ui.T->item(i,0)->text().toInt()] = bc;
  };
};

