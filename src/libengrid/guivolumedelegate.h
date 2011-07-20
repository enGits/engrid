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
#ifndef GUIVOLUMEDELEGATE_H
#define GUIVOLUMEDELEGATE_H

#include <QItemDelegate>

class GuiVolumeDelegate : public QItemDelegate
{

  Q_OBJECT;

public:

  void     paint        (QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;
  QWidget* createEditor (QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
  void     setEditorData(QWidget *editor, const QModelIndex &index) const;
  void     setModelData (QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
  void     setFirstCol  (int c) { first_column = c; }

private slots:

  void commitAndCloseEditor();

private:

  int first_column;

};

#endif // GUIVOLUMEDELEGATE_H
