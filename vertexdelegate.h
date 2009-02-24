#ifndef VERTEXDELEGATE_H
#define VERTEXDELEGATE_H

#include <QItemDelegate>

class VertexDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    VertexDelegate(int durationColumn, QList<QString> list, QObject *parent = 0);

    void paint(QPainter *painter, const QStyleOptionViewItem &option,
               const QModelIndex &index) const;
    QWidget *createEditor(QWidget *parent,
                          const QStyleOptionViewItem &option,
                          const QModelIndex &index) const;
    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const;

private slots:
    void commitAndCloseEditor();

private:
	QList<QString> list;
    int durationColumn;
};

#endif
