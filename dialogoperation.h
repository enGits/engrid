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
#ifndef dialogoperation_H
#define dialogoperation_H

template <class UI>
class DialogOperation;

#include "operation.h"

#include <QDialog>
#include <QListWidget>
#include <QTextStream>

#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>

template <class UI>
class DialogOperation : public QDialog, 
                        public Operation
{
  
protected: // attributes
  
  /** The user interface definition from QtDesigner */
  UI ui;
  
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
   * Get a set qith all seleceted items from a QListWidget.
   * @param lw  The QListWidget.
   * @param sel On return, this will hold all items.
   */
  void getSelectedItems(QListWidget *lw, QSet<int> &sel);
  
  /**
   * Get a set qith all seleceted items from a QListWidget.
   * @param lw  The QListWidget.
   * @param sel On return, this will hold all items.
   */
  void getSelectedItems(QListWidget *lw, QSet<QString> &sel);
  
  /**
   * Fill a QListWidget with all available boundary codes from a grid.
   * @param lw   The QListWidget to fill.
   * @param grid The grid to use.
   */
  void populateBoundaryCodes(QListWidget *lw, vtkUnstructuredGrid *grid);
  
  virtual void before() {};
  virtual void operator()();
  
  //connect(const QObject* a, const char* b, const QObject* c, const char* d) { QObject::connect(a,b,c,d); };
};


template <class UI>
DialogOperation<UI>::DialogOperation()
{
  ui.setupUi(this);
};

template <class UI>
template <class T>
void DialogOperation<UI>::addListItem(QListWidget *lw, T item, bool checked)
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

template <class UI>
template <class T>
bool DialogOperation<UI>::checkListItem(QListWidget *lw, T item)
{
  QString text = "";
  QTextStream ts(&text);
  ts << item;
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->text() == text) {
      if (lw->item(i)->checkState() == Qt::Checked) return true;
    };
  };
  return false;
};

template <class UI>
void DialogOperation<UI>::getSelectedItems(QListWidget *lw, QSet<QString> &sel)
{
  sel.clear();
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->checkState() == Qt::Checked) {
      QString item = lw->item(i)->text();
      sel.insert(item);
    };
  };
};

template <class UI>
void DialogOperation<UI>::getSelectedItems(QListWidget *lw, QSet<int> &sel)
{
  sel.clear();
  for (int i = 0; i < lw->count(); ++i) {
    if (lw->item(i)->checkState() == Qt::Checked) {
      int item = lw->item(i)->text().toInt();
      sel.insert(item);
    };
  };
};

template <class UI>
void DialogOperation<UI>::populateBoundaryCodes(QListWidget *lw, vtkUnstructuredGrid *grid)
{
  try {
    QSet<int> bcodes;
    EG_VTKDCC(vtkIntArray,cell_code,grid,"cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      int ct = grid->GetCellType(i);
      if ((ct == VTK_TRIANGLE) || (ct == VTK_QUAD)) {
        bcodes.insert(cell_code->GetValue(i));
      };
    };
    int bc;
    foreach(bc, bcodes) addListItem(lw, bc, false);
  } catch (Error err) {
    err.display();
  };
};

template <class UI>
void DialogOperation<UI>::operator()()
{
  bool ok = true;
  try {
    checkGrid();
    before();
  } catch (Error err) {
    err.display();
    ok = false;
  };
  if (ok) {
    if (QDialog::exec()) {
      try {
        Operation::operator()();
      } catch (Error err) {
        err.display();
      };
    };
  };
};

#endif
