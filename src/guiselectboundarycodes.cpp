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
#include "guiselectboundarycodes.h"
#include <vtkIntArray.h>
#include <vtkCellData.h>

#include "guimainwindow.h"

GuiSelectBoundaryCodes::GuiSelectBoundaryCodes()
{
  disableAutoSet();
}

void GuiSelectBoundaryCodes::setDisplayBoundaryCodes(const QSet<int> &bcs)
{
  m_DisplayBoundaryCodes = bcs;
}

void GuiSelectBoundaryCodes::before()
{
  populateBoundaryCodes(ui.listWidget);
  int row = 0;
  for (QSet<int>::iterator i = m_BoundaryCodes.begin(); i != m_BoundaryCodes.end(); ++i) {
    bool checked = m_DisplayBoundaryCodes.contains(*i);
    if (checked) {
      ui.listWidget->item(row)->setCheckState(Qt::Checked);
    } else {
      ui.listWidget->item(row)->setCheckState(Qt::Unchecked);
    }
    ++row;
  }
  connect(ui.pushButtonSelect, SIGNAL(clicked()), this, SLOT(selectAll()));
  connect(ui.pushButtonDeselect, SIGNAL(clicked()), this, SLOT(deselectAll()));
  connect(ui.pushButton_SaveSelectionAsGrid, SIGNAL(clicked()), this, SLOT(saveSelectionAsGrid()));
}

void GuiSelectBoundaryCodes::operate()
{
  getSelectedItems(ui.listWidget, m_DisplayBoundaryCodes);
}

void GuiSelectBoundaryCodes::selectAll()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiSelectBoundaryCodes::deselectAll()
{
  for (int i = 0; i < ui.listWidget->count(); ++i) {
    ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  }
}

void GuiSelectBoundaryCodes::saveSelectionAsGrid()
{
  getSelectedItems(ui.listWidget, m_DisplayBoundaryCodes);
  QVector <vtkIdType> selected_cells;
  getSurfaceCells(m_DisplayBoundaryCodes, selected_cells, m_Grid);
  writeCells(m_Grid, selected_cells, GuiMainWindow::pointer()->getFilePath() + "selection.vtu" );
}

void GuiSelectBoundaryCodes::getSelectedBoundaryCodes(QSet<int> &bcs)
{
  bcs = m_DisplayBoundaryCodes;
}


