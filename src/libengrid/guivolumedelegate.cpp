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
#include "guivolumedelegate.h"

#include <QComboBox>
#include <QtDebug>

void GuiVolumeDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
  if (index.column() >= first_column) {
    QString text = index.model()->data(index, Qt::DisplayRole).toString();
    QStyleOptionViewItem myOption = option;
    myOption.displayAlignment = Qt::AlignRight | Qt::AlignVCenter;
    drawDisplay(painter, myOption, myOption.rect, text);
    drawFocus(painter, myOption, myOption.rect);
  } else{
    QItemDelegate::paint(painter, option, index);
  }
}

QWidget *GuiVolumeDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
  if (index.column() >= first_column) {
    QComboBox *comboBox = new QComboBox(parent);
    comboBox->addItem("green");
    comboBox->addItem("yellow");
    comboBox->addItem(" ");
    connect(comboBox, SIGNAL(currentIndexChanged ( int )), this, SLOT(commitAndCloseEditor()));
    return comboBox;
  } else {
    return QItemDelegate::createEditor(parent, option, index);
  }
}

void GuiVolumeDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const
{
  if (index.column() >= first_column) {
    int secs = index.model()->data(index, Qt::DisplayRole).toInt();
    QComboBox *comboBox = qobject_cast<QComboBox *>(editor);
    comboBox->setCurrentIndex(secs);
  } else {
    QItemDelegate::setEditorData(editor, index);
  }
}

void GuiVolumeDelegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
  if (index.column() >= first_column) {
    QComboBox *comboBox = qobject_cast<QComboBox *>(editor);
    model->setData(index, comboBox->currentText());
  } else {
    QItemDelegate::setModelData(editor, model, index);
  }
}

void GuiVolumeDelegate::commitAndCloseEditor()
{
  qDebug()<<"commitAndCloseEditor called";
  QComboBox *editor = qobject_cast<QComboBox *>(sender());
  emit commitData(editor);
  emit closeEditor(editor);
}

