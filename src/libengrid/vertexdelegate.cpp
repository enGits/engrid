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
#include <QtGui>
#include <QComboBox>

#include "vertexdelegate.h"
#include "surfacemesher.h"

#include <iostream>
using namespace std;

VertexDelegate::VertexDelegate(int Column, QList<QString> list, QObject *parent)
    : QItemDelegate(parent)
{
    this->list = list;
    this->Column = Column;
}


void VertexDelegate::paint(QPainter *painter,
                          const QStyleOptionViewItem &option,
                          const QModelIndex &index) const
{
    if (index.column() == Column) {
        QString text = index.model()->data(index, Qt::DisplayRole).toString();

        QStyleOptionViewItem myOption = option;
        myOption.displayAlignment = Qt::AlignRight | Qt::AlignVCenter;

        drawDisplay(painter, myOption, myOption.rect, text);
        drawFocus(painter, myOption, myOption.rect);
    } else{
        QItemDelegate::paint(painter, option, index);
    }
}

QWidget *VertexDelegate::createEditor(QWidget *parent,
        const QStyleOptionViewItem &option,
        const QModelIndex &index) const
{
    if (index.column() == Column) {
        QComboBox *ComboEdit = new QComboBox(parent);
	foreach(QString str,list) ComboEdit->addItem(str);
        connect(ComboEdit, SIGNAL(currentIndexChanged ( int )),this, SLOT(commitAndCloseEditor()));
        return ComboEdit;
    } else {
        return QItemDelegate::createEditor(parent, option, index);
    }
}

void VertexDelegate::setEditorData(QWidget *editor,
                                  const QModelIndex &index) const
{
    if (index.column() == Column) {
//         int secs = index.model()->data(index, Qt::DisplayRole).toInt();
//         QComboBox *ComboEdit = qobject_cast<QComboBox *>(editor);
// 	ComboEdit->setCurrentIndex(secs);
    } else {
        QItemDelegate::setEditorData(editor, index);
    }
}

void VertexDelegate::setModelData(QWidget *editor,
                                 QAbstractItemModel *model,
                                 const QModelIndex &index) const
{
    if (index.column() == Column) {
        QComboBox *ComboEdit = qobject_cast<QComboBox *>(editor);
        model->setData(index, ComboEdit->currentText());
    } else {
        QItemDelegate::setModelData(editor, model, index);
    }
}

void VertexDelegate::commitAndCloseEditor()
{
    QComboBox *editor = qobject_cast<QComboBox *>(sender());
    emit commitData(editor);
    emit closeEditor(editor);
}
