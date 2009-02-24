#include <QtGui>

#include "vertexdelegate.h"

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
        QComboBox *timeEdit = new QComboBox(parent);
	foreach(QString str,list) timeEdit->addItem(str);
        connect(timeEdit, SIGNAL(editingFinished()),
                this, SLOT(commitAndCloseEditor()));
        return timeEdit;
    } else {
        return QItemDelegate::createEditor(parent, option, index);
    }
}

void VertexDelegate::setEditorData(QWidget *editor,
                                  const QModelIndex &index) const
{
	cout<<"DDD"<<endl;

    if (index.column() == durationColumn) {
        int secs = index.model()->data(index, Qt::DisplayRole).toInt();
        QComboBox *timeEdit = qobject_cast<QComboBox *>(editor);
	timeEdit->setCurrentIndex(secs);
    } else {
        QItemDelegate::setEditorData(editor, index);
    }
}

void VertexDelegate::setModelData(QWidget *editor,
                                 QAbstractItemModel *model,
                                 const QModelIndex &index) const
{
    if (index.column() == durationColumn) {
        QComboBox *timeEdit = qobject_cast<QComboBox *>(editor);
        model->setData(index, timeEdit->currentText());
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
