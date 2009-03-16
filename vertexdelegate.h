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
