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
#ifndef operation_H
#define operation_H

class Operation;
class GuiMainWindow;

#include "egvtkobject.h"
#include "vertexmeshdensity.h"
#include "meshpartition.h"

#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>

#include <QThread>
#include <QMutex>
#include <QListWidget>

#include <typeinfo>

class OperationThread : public QThread
{
  
private:
  
  Operation *op;

  static QVector<vtkIdType>     m_static_DummyCells;  /// dummy to initialise references
  static QVector<int>           m_static_DummyRCells; /// dummy to initialise references
  static QVector<vtkIdType>     m_static_DummyNodes;  /// dummy to initialise references
  static QVector<int>           m_static_DummyRNodes; /// dummy to initialise references
  static QVector<QVector<int> > m_static_DummyN2N;    /// dummy to initialise references
  static QVector<QVector<int> > m_static_DummyN2C;    /// dummy to initialise references
  static QVector<QVector<int> > m_static_DummyC2C;    /// dummy to initialise references

protected:
  
  virtual void run();
  
public:
  
  void setOperation(Operation *an_op) { op = an_op; }
  
};

/**
 * This is the base class for all mesh operations.
 * Operations will typically be triggered by a Qt event; the MainWindow
 * object will call operator() with the current grid as parameter.
 */
class Operation : public EgVtkObject
{
  
  friend class OperationThread;
  OperationThread thread;
  
private: // static attributes
  
  static QSet<Operation*> garbage_operations;
  
private: // attributes
  
  bool               gui;
  bool               m_quicksave;             /// save grid after operation finished?
  bool               m_resetoperationcounter; /// reset operation counter after operation finished? (default is false)
  bool               autoset;
  Error             *err;
  QString            volume_name;
  
private: // methods
  
  //void initMapping();
  
protected: // attributes
  
  vtkUnstructuredGrid    *grid;   /// The main grid the operation operates on.
  MeshPartition           m_Part; /// the partition containing the subset of cells and nodes
  QVector<vtkIdType>     &cells;  /// reference for backwards compatibility
  QVector<int>           &_cells; /// reference for backwards compatibility
  QVector<vtkIdType>     &nodes;  /// reference for backwards compatibility
  QVector<int>           &_nodes; /// reference for backwards compatibility
  QVector<QVector<int> > &n2n;    /// reference for backwards compatibility
  QVector<QVector<int> > &n2c;    /// reference for backwards compatibility
  QVector<QVector<int> > &c2c;    /// reference for backwards compatibility

protected: // methods
  
  void checkGrid();
  void updateActors();
  GuiMainWindow* mainWindow();
  virtual void operate() = 0;

public: // methods
  
  Operation();
  virtual ~Operation();
  void del();
    
  void setGrid(vtkUnstructuredGrid *ug) { grid = ug; }
  void setAllCells();
  void setAllVolumeCells();
  void setAllSurfaceCells();
  void setVolume(QString volume_name);
  template <class T> void setCells(const T &cls);
  void setMeshPartition(const MeshPartition &part);

  void setGui() { gui = true; }
  OperationThread& getThread() { return thread; }
  void enableAutoSet() { autoset = true; }
  void disableAutoSet() { autoset = false; }
  void setQuickSave(bool b) { m_quicksave = b; }
  void setResetOperationCounter(bool b) { m_resetoperationcounter=b; }
  
  /**
   * Fill a QListWidget with all available boundary codes from a grid.
   * @param lw   The QListWidget to fill.
   * @param grid The grid to use.
   */
  void populateBoundaryCodes(QListWidget *lw);
  
  /**
   * Fill a QListWidget with all available volumes from a grid.
   * @param lw   The QListWidget to fill.
   * @param grid The grid to use.
   */
  void populateVolumes(QListWidget *lw);

  virtual void operator()();

  static void collectGarbage();
};

//End of class Operation




template <class T>
void Operation::setCells(const T &cls)
{
  m_Part.setGrid(grid);
  m_Part.setCells(cls);
  cells = m_Part.getCells();
  nodes = m_Part.getNodes();
  _cells = m_Part.getLocalCells();
  _nodes = m_Part.getLocalNodes();
  n2n = m_Part.getN2N();
  n2c = m_Part.getN2C();
  c2c = m_Part.getC2C();
}

#endif
