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
  populateBoundaryCodes(m_Ui.listWidget);
  for (int row = 0; row < m_Ui.listWidget->count(); ++row) {
    QString txt = m_Ui.listWidget->item(row)->text();
    QStringList words = txt.split(":");
    if (words.size() < 1) {
      EG_BUG;
    }
    int bc = words[0].toInt();
    bool checked = m_DisplayBoundaryCodes.contains(bc);
    if (checked) {
      m_Ui.listWidget->item(row)->setCheckState(Qt::Checked);
    } else {
      m_Ui.listWidget->item(row)->setCheckState(Qt::Unchecked);
    }
  }

  connect(m_Ui.pushButtonSelect, SIGNAL(clicked()), this, SLOT(selectAll()));
  connect(m_Ui.pushButtonDeselect, SIGNAL(clicked()), this, SLOT(deselectAll()));
  connect(m_Ui.pushButton_SaveSelectionAsGrid, SIGNAL(clicked()), this, SLOT(saveSelectionAsGrid()));
}

void GuiSelectBoundaryCodes::operate()
{
  getSelectedItems(m_Ui.listWidget, m_DisplayBoundaryCodes);
}

void GuiSelectBoundaryCodes::selectAll()
{
  for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
    m_Ui.listWidget->item(i)->setCheckState(Qt::Checked);
  }
}

void GuiSelectBoundaryCodes::deselectAll()
{
  for (int i = 0; i < m_Ui.listWidget->count(); ++i) {
    m_Ui.listWidget->item(i)->setCheckState(Qt::Unchecked);
  }
}

void GuiSelectBoundaryCodes::saveSelectionAsGrid()
{
  getSelectedItems(m_Ui.listWidget, m_DisplayBoundaryCodes);
  QVector <vtkIdType> selected_cells;
  getSurfaceCells(m_DisplayBoundaryCodes, selected_cells, m_Grid);
  writeCells(m_Grid, selected_cells, GuiMainWindow::pointer()->getFilePath() + "selection.vtu" );
}

void GuiSelectBoundaryCodes::getSelectedBoundaryCodes(QSet<int> &bcs)
{
  bcs = m_DisplayBoundaryCodes;
}


