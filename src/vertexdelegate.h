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
#ifndef VERTEXDELEGATE_H
#define VERTEXDELEGATE_H

#include <QItemDelegate>
#include <QTableWidgetItem>

//Thanks to jpnurmi from qt@irc.freenode.net for this class ;)
class TriStateTableWidgetItem : public QTableWidgetItem
{
public:
  TriStateTableWidgetItem() : QTableWidgetItem()
  {
    QTableWidgetItem::setCheckState(Qt::Checked);
    setData(Qt::CheckStateRole,Qt::Checked);
  }
  void setData(int role, const QVariant& value)
  {
    QVariant tmp = value;
    if (role == Qt::CheckStateRole)
    {
      if (data(role) == Qt::Unchecked && value == Qt::Checked)
        tmp = Qt::PartiallyChecked;
      else if (data(role) == Qt::PartiallyChecked && value == Qt::Checked)
        tmp = Qt::Checked;
      else
        tmp = Qt::Unchecked;
    }
    QTableWidgetItem::setData(role, tmp);
  }
  void setCheckState(Qt::CheckState S){
    QTableWidgetItem::setData(Qt::CheckStateRole, S);
  }
};

/// Creates a combobox for use in a table
class VertexDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    VertexDelegate(int Column, QList<QString> list, QObject *parent = 0);

    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;
    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;

private slots:
    void commitAndCloseEditor();

private:
    QList<QString> list;
    int Column;
};

#endif
