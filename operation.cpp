//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#include "operation.h"
#include "guimainwindow.h"

#include <QApplication>


QSet<Operation*> Operation::garbage_operations;

void Operation::collectGarbage()
{
  QSet<Operation*> delete_operations;
  foreach (Operation *op, garbage_operations) {
    if (!op->getThread().isRunning()) {
      delete_operations.insert(op);
      cout << "deleting Operation " << op << endl;
      delete op;
    };
  };
  foreach (Operation *op, delete_operations) {
    garbage_operations.remove(op);
  };
};

Operation::Operation()
{
  grid = NULL;
  gui = false;
  err = NULL;
  autoset = true;
};

Operation::~Operation()
{
  if (err) {
    err->display();
    delete err;
  };
};

void OperationThread::run()
{
  try {
    GuiMainWindow::lock();
    GuiMainWindow::pointer()->setBusy();
    op->operate();
  } catch (Error err) {
    op->err = new Error();
    *(op->err) = err;
  };
  GuiMainWindow::unlock();
  GuiMainWindow::pointer()->setIdle();
};

void Operation::operator()()
{
  if (gui) {
    if (GuiMainWindow::tryLock()) {
      checkGrid();
      thread.setOperation(this);
      GuiMainWindow::unlock();
      thread.start(QThread::LowPriority);
    } else {
      QMessageBox::warning(NULL, "not permitted", "Operation is not permitted while background process is running!");
    };
  } else {
    checkGrid();
    operate();
  };
  /*
  checkGrid();
  operate();
  */
};

void Operation::setAllCells()
{
  QVector<vtkIdType> all_cells;
  getAllCells(all_cells, grid);
  setCells(all_cells);
};

void Operation::setAllVolumeCells()
{
  QVector<vtkIdType> cells;
  getAllVolumeCells(cells, grid);
  setCells(cells);
};

void Operation::setAllSurfaceCells()
{
  QVector<vtkIdType> cells;
  getAllSurfaceCells(cells, grid);
  setCells(cells);
};

void Operation::initMapping()
{
  nodes_map.resize(nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    nodes_map[i_nodes] = nodes[i_nodes];
  };
  cells_map.resize(cells.size());
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    cells_map[i_cells] = cells[i_cells];
  };
};

void Operation::checkGrid()
{
  if (grid == NULL) {
    grid = GuiMainWindow::pointer()->getGrid();
  };
  if ((cells.size() == 0) && autoset) {
    setAllCells();
  };
};

void Operation::updateActors()
{
  mainWindow()->updateActors();
};

GuiMainWindow* Operation::mainWindow()
{
  return GuiMainWindow::pointer();
};
