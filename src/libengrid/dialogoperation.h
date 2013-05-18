// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#ifndef dialogoperation_H
#define dialogoperation_H

template <class UI, class OP>
class DialogOperation;

#include "operation.h"
// #include "guimainwindow.h"

#include <QDialog>
#include <QListWidget>
#include <QTextStream>

#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

template <class UI, class OP>
class DialogOperation : public QDialog, 
                        public OP
{
  
protected: // attributes
  
  /** The user interface definition from QtDesigner */
  UI m_Ui;

protected: // methods

  EgVtkObject::l2g_t getPartNodes()      { return this->m_Part.getNodes(); }
  EgVtkObject::l2g_t getPartCells()const { return this->m_Part.getCells(); }
  EgVtkObject::g2l_t getPartLocalNodes() { return this->m_Part.getLocalNodes(); }
  EgVtkObject::g2l_t getPartLocalCells() { return this->m_Part.getLocalCells(); }
  EgVtkObject::l2l_t getPartN2N()        { return this->m_Part.getN2N(); }
  EgVtkObject::l2l_t getPartN2C()        { return this->m_Part.getN2C(); }
  EgVtkObject::l2l_t getPartC2C()        { return this->m_Part.getC2C(); }
  
public: // methods
  
  /** Generic constructor to set up the user interface */
  DialogOperation();
  
  /**
   * Create a checkable QListWidgetItem and add it to a QListWidget.
   * @param lw the QListWidget to add the item to
   * @param item the text of the item
   * @param checked the status of the item checked/unchecked
   */
  template <class T>
  void addListItem(QListWidget *lw, T item, bool checked = false);
  
  /**
   * Check if a certain text item is checked within a QListWidget.
   * @param lw the QListWidget to check
   * @param item the item to check
   * @return true if the item is checked and false if not
   */
  template <class T>
  bool checkListItem(QListWidget *lw, T item);

  /**
   * Get the name of the selected volume.
   * @param lw the QListWidget for volume selection
   * @return the name of the selected volume
   */
  QString getSelectedVolume(QListWidget *lw);

  /**
   * Get a set with all seleceted items from a QListWidget.
   * @param lw  The QListWidget.
   * @param sel On return, this will hold all items.
   */
  void getSelectedItems(QListWidget *lw, QSet<int> &sel);
  
  /**
   * Get a set with all seleceted items from a QListWidget.
   * @param lw  The QListWidget.
   * @param sel On return, this will hold all items.
   */
  void getSelectedItems(QListWidget *lw, QSet<QString> &sel);
  
  virtual void before() {}
  virtual void operator()();
  
  //connect(const QObject* a, const char* b, const QObject* c, const char* d) { QObject::connect(a,b,c,d); };
};


template <class UI, class OP>
DialogOperation<UI,OP>::DialogOperation()
{
  m_Ui.setupUi(this);
};

template <class UI, class OP>
template <class T>
void DialogOperation<UI,OP>::addListItem(QListWidget *lw, T item, bool checked)
{
  QListWidgetItem *lwi = new QListWidgetItem(lw);
  if (checked) lwi->setCheckState(Qt::Checked);
  else         lwi->setCheckState(Qt::Unchecked);
  QString text = "";
  QTextStream ts(&text);
  ts << item;
  lwi->setText(text);
  lwi->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
};

template <class UI, class OP>
template <class T>
bool DialogOperation<UI,OP>::checkListItem(QListWidget *lw, T item)
{
  QString text = "";
  QTextStream ts(&text);
  ts << item;
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->text() == text) {
      if (lw->item(i)->checkState() == Qt::Checked) return true;
    }
  }
  return false;
}

template <class UI, class OP>
QString DialogOperation<UI,OP>::getSelectedVolume(QListWidget *lw)
{
  QString volume_name;
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->isSelected()) {
      volume_name = lw->item(i)->text();
    }
  }
  return volume_name;
}

template <class UI, class OP>
void DialogOperation<UI,OP>::getSelectedItems(QListWidget *lw, QSet<QString> &sel)
{
  sel.clear();
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->checkState() == Qt::Checked) {
      QString item = lw->item(i)->text();
      sel.insert(item);
    }
  }
}

template <class UI, class OP>
void DialogOperation<UI,OP>::getSelectedItems(QListWidget *lw, QSet<int> &sel)
{
  sel.clear();
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->checkState() == Qt::Checked) {
      QString item_txt = lw->item(i)->text();
      QStringList items = item_txt.split(":");
      int item = items[0].toInt(); ///\todo UiUiUi
      sel.insert(item);
    }
  }
}

template <class UI, class OP>
void DialogOperation<UI,OP>::operator()()
{
  bool ok = true;
  //Prepare the GUI
  try {
    OP::checkGrid();
    before();
  } catch (Error err) {
    err.display();
    ok = false;
  };
  //Run the GUI
  if (ok) {
    try {
      if (!QDialog::exec()) {
        ok = false;
      }
    } catch (Error err) {
      err.display();
    }
  }
  //Run the operation
  if (ok) {
    try {
      OP::operator()();
    } catch (Error err) {
      err.display();
    }
  }
}

#endif
