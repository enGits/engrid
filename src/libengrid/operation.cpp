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
#include "operation.h"

#include "guimainwindow.h"
#include "egvtkobject.h"
#include "setboundarycode.h"

#include <vtkTriangleFilter.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkCharArray.h>

#include <QApplication>
#include <QTime>

#include "geometrytools.h"

using namespace GeometryTools;

QSet<Operation*> Operation::garbage_operations;

QVector<vtkIdType>     m_static_DummyCells;
QVector<int>           m_static_DummyRCells;
QVector<vtkIdType>     m_static_DummyNodes;
QVector<int>           m_static_DummyRNodes;
QVector<QVector<int> > m_static_DummyN2N;
QVector<QVector<int> > m_static_DummyN2C;
QVector<QVector<int> > m_static_DummyC2C;

void Operation::collectGarbage()
{
  QSet<Operation*> delete_operations;

  foreach (Operation *op, garbage_operations)
  {
    if (!op->getThread().isRunning()) {
      delete_operations.insert(op);
      cout << "deleting Operation " << op << endl;
      delete op;
    }
  }

  foreach (Operation *op, delete_operations)
  {
    garbage_operations.remove(op);
  }
}

Operation::Operation()
{
  m_Grid = NULL;
  lock_gui = false;
  m_quicksave = false;
  m_resetoperationcounter = false;
  err = NULL;
  autoset = true;
  m_TypeName = "undefined";
  m_MenuText = "undefined";
  m_Verbose = false;
}

Operation::~Operation()
{
  if (err) {
    err->display();
    delete err;
  }
}

void Operation::del() 
{ 
  garbage_operations.insert(this); 
}

void OperationThread::run()
{
  try {
    GuiMainWindow::lock();
    GuiMainWindow::pointer()->setBusy();
    op->operate();
    cout << "secs. for " << qPrintable(op->getTypeName()) << ": " << op->elapsedTime() << endl;
  } catch (Error err) {
    op->err = new Error();
    *(op->err) = err;
  }
  GuiMainWindow::unlock();
  GuiMainWindow::pointer()->setIdle();
}

void Operation::setTypeName(QString name)
{
  int i = 0;
  while ((i < name.size()) && (name[i].isDigit())) {
    ++i;
  }
  m_TypeName = name.right(name.size() - i);
}

void Operation::operator()()
{
  setStartTime();
  if (lock_gui) {
    if (GuiMainWindow::tryLock()) {
      checkGrid();
      thread.setOperation(this);
      GuiMainWindow::unlock();
      thread.start(QThread::LowPriority);
    } else {
      QMessageBox::warning(NULL, "not permitted", "Operation is not permitted while background process is running!");
    }
  } else {
    checkGrid();
    const bool gui_thread = QThread::currentThread() == QCoreApplication::instance()->thread();
    if (gui_thread) {
      try {
        QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
        operate();
        QApplication::restoreOverrideCursor();
        cout << "secs. for " << qPrintable(getTypeName()) << ": " << elapsedTime() << endl;
      } catch (Error err) {
        err.display();
      }
    } else {
      operate();
    }
    if(m_resetoperationcounter) GuiMainWindow::pointer()->resetOperationCounter();
    if(m_quicksave) GuiMainWindow::pointer()->quickSave();
  }
}

void Operation::setAllCells()
{
  QVector<vtkIdType> all_cells;
  getAllCells(all_cells, m_Grid);
  setCells(all_cells);
}

void Operation::setAllVolumeCells()
{
  QVector<vtkIdType> cells;
  getAllVolumeCells(cells, m_Grid);
  setCells(cells);
}

void Operation::setAllSurfaceCells()
{
  QVector<vtkIdType> cells;
  getAllSurfaceCells(cells, m_Grid);
  setCells(cells);
}

void Operation::setVolume(QString volume_name)
{
  m_Part.setGrid(m_Grid);
  m_Part.setVolume(volume_name);
}

void Operation::setMeshPartition(const MeshPartition &part)
{
  m_Part.setGrid(part.getGrid());
  m_Part.setCells(part.getCells());
  m_Grid = m_Part.getGrid();
}

void Operation::checkGrid()
{
  if (m_Grid == NULL) {
    m_Grid = GuiMainWindow::pointer()->getGrid();
  }
  l2g_t cells = getPartCells();
  if ((cells.size() == 0) && autoset) {
    setAllCells();
  }
}

void Operation::updateActors()
{
  mainWindow()->updateActors();
}

GuiMainWindow* Operation::mainWindow()
{
  return GuiMainWindow::pointer();
}

void Operation::populateBoundaryCodes(QListWidget *lw)
{
  QVector<int> bcs;
  mainWindow()->getAllBoundaryCodes(bcs);
  foreach(int bc, bcs) {
    QListWidgetItem *lwi = new QListWidgetItem(lw);
    lwi->setCheckState(Qt::Unchecked);
    QString text = "";
    QTextStream ts(&text);
    ts << bc << " : " << GuiMainWindow::pointer()->getBC(bc).getName();
    lwi->setText(text);
    lwi->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  }
}

void Operation::populateVolumes(QListWidget *lw)
{
  QList<VolumeDefinition> vols = mainWindow()->getAllVols();
  foreach (VolumeDefinition V, vols) {
    QListWidgetItem *lwi = new QListWidgetItem(lw);
    lwi->setText(V.getName());
  }
}

void Operation::eliminateDuplicateCells(bool surf_only)
{
  QVector<QVector<vtkIdType> > cell_nodes(m_Grid->GetNumberOfCells());

  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (!surf_only || isSurface(id_cell, m_Grid)) {
      EG_GET_CELL(id_cell, m_Grid);
      QVector<vtkIdType> nodes(num_pts);
      for (int i = 0; i < num_pts; ++i) {
        nodes[i] = pts[i];
      }
      qSort(nodes);
      cell_nodes[id_cell] = nodes;
    }
  }
  QList<vtkIdType> new_cells;
  for (vtkIdType id_cell1 = 0; id_cell1 < m_Grid->GetNumberOfCells(); ++id_cell1) {
    bool duplicate_cell = false;
    if (!surf_only || isSurface(id_cell1, m_Grid)) {
      EG_GET_CELL(id_cell1, m_Grid);
      for (int i = 0; i < num_pts; ++i) {
        for (int j = 0; j < m_Part.n2cGSize(pts[i]); ++j) {
          vtkIdType id_cell2 = m_Part.n2cGG(pts[i], j);
          if (id_cell1 != id_cell2) {
            if (!surf_only || isSurface(id_cell2, m_Grid)) {
              if (cell_nodes[id_cell1] == cell_nodes[id_cell2]) {
                duplicate_cell = true;
                break;
              }
            }
          }
        }
      }
    }
    if (!duplicate_cell) {
      new_cells.append(id_cell1);
    }
  }
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  makeCopy(m_Grid, new_grid, new_cells);
  makeCopy(new_grid, m_Grid);
}

void Operation::createFeatureBcs(double feature_angle)
{
  m_New2OldBc.clear();
  SetBoundaryCode set_bc;
  set_bc.setGrid(m_Grid);
  set_bc.setAllSurfaceCells();
  QSet <int> all_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  foreach (int bc, all_bcs) {
    m_New2OldBc[bc] = bc;
  }

  QMap<int, QList<vtkIdType> > oldbc2cells;
  int max_bc = 1;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    int bc = cell_code->GetValue(id_cell);
    max_bc = max(max_bc, bc + 1);
    if (all_bcs.contains(bc)) {
      oldbc2cells[bc] << id_cell;
    }
  }
  bool done = false;
  int old_max_bc = max_bc;
  do {
    vtkIdType id_start = -1;
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      if (cell_code->GetValue(id_cell) < old_max_bc && isSurface(id_cell, m_Grid)) {
        id_start = id_cell;
        set_bc.setOLdBC(cell_code->GetValue(id_cell));
        break;
      }
    }
    if (id_start == -1) {
      done = true;
    } else {
      set_bc.setFeatureAngle(feature_angle);
      set_bc.setNewBC(max_bc);
      set_bc.setProcessAll(true);
      set_bc.setSelectAllVisible(false);
      set_bc.setOnlyPickedCell(false);
      set_bc.setOnlyPickedCellAndNeighbours(false);
      set_bc.setStart(id_start);
      set_bc();
      ++max_bc;
    }
  } while (!done);
  m_New2OldBc.clear();
  foreach (int bc, oldbc2cells.keys()) {
    foreach (vtkIdType id_cell, oldbc2cells[bc]) {
      m_New2OldBc[cell_code->GetValue(id_cell)] = bc;
    }
  }
}

void Operation::restoreNormalBcs()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    int bc = cell_code->GetValue(id_cell);
    cell_code->SetValue(id_cell, m_New2OldBc[bc]);
  }
}

