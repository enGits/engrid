#include <QtGui>

#include "vertexdelegate.h"
#include "createspecialmapping.h"

#include <iostream>
using namespace std;

VertexDelegate::VertexDelegate(int durationColumn, QList<QString> list, QObject *parent)
    : QItemDelegate(parent)
{
    this->list = list;
    this->durationColumn = durationColumn;
}

void VertexDelegate::paint(QPainter *painter,
                          const QStyleOptionViewItem &option,
                          const QModelIndex &index) const
{
    if (index.column() == durationColumn) {
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
    if (index.column() == durationColumn) {
        QComboBox *ComboEdit = new QComboBox(parent);
	foreach(QString str,list) ComboEdit->addItem(str);
        connect(ComboEdit, SIGNAL(editingFinished()),
                this, SLOT(commitAndCloseEditor()));
        return ComboEdit;
    } else {
        return QItemDelegate::createEditor(parent, option, index);
    }
}

void VertexDelegate::setEditorData(QWidget *editor,
                                  const QModelIndex &index) const
{
    if (index.column() == durationColumn) {
//         int secs = index.model()->data(index, Qt::DisplayRole).toInt();
        QComboBox *ComboEdit = qobject_cast<QComboBox *>(editor);
// 	ComboEdit->setCurrentIndex(secs);
    } else {
        QItemDelegate::setEditorData(editor, index);
    }
}

void VertexDelegate::setModelData(QWidget *editor,
                                 QAbstractItemModel *model,
                                 const QModelIndex &index) const
{
    if (index.column() == durationColumn) {
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
