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
#ifndef operation_H
#define operation_H

class Operation;
class GuiMainWindow;

#include "egvtkobject.h"
#include "vertexmeshdensity.h"
#include "meshpartition.h"
#include "timer.h"

#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>

#include <QThread>
#include <QMutex>
#include <QListWidget>
#include <QTime>

#include <typeinfo>

#define EG_TYPENAME setTypeName(QString(typeid(*this).name()))

class OperationThread : public QThread
{
  
private:
  
  Operation *op;

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
  
  /// If true, attempts to lock the mutex when Operation::operator()() is called. If the lock was obtained, all other "lock_gui operations" will not work until the current "lock_gui opration" is done.
  bool               lock_gui;
  
  bool               m_quicksave;             ///< save grid after operation finished?
  bool               m_resetoperationcounter; ///< reset operation counter after operation finished? (default is false)
  bool               autoset;
  Error             *err;
  QString            volume_name;
  QString            m_TypeName;              ///< typename retrieved from typeid(this).name()
  QTime              m_StartTime;             ///< start time for run-time information

private: // methods

  void setStartTime() { m_StartTime = QTime::currentTime(); }
  int  elapsedTime() { return m_StartTime.secsTo(QTime::currentTime()); }
  
protected: // attributes
  
  vtkUnstructuredGrid* m_Grid;     ///< The main grid the operation operates on.
  vtkUnstructuredGrid* m_RestGrid; ///< The remainder grid (not part of the selected volume)
  MeshPartition        m_Part;     ///< the partition containing the subset of cells and nodes
  Timer                m_Timer;    ///< Timer object for periodic output
  QString              m_MenuText; ///< The menu entry (mainly for plugins)
  bool                 m_Verbose;  ///< General flag to trigger more output

protected: // methods
  
  void checkGrid();
  void updateActors();
  GuiMainWindow* mainWindow();
  virtual void operate() = 0;
  void setTypeName(QString name);

  /**
   * Eliminate cells with identical node indices.
   * @param surf_only if set to false all cells will be checked, otherwise surface cells only.
   *                  Careful with eliminating volume cells -- this can be extremely slow.
   */
  void eliminateDuplicateCells(bool surf_only = true);

  l2g_t getPartNodes()       { return m_Part.getNodes(); }
  l2g_t getPartCells() const { return m_Part.getCells(); }
  g2l_t getPartLocalNodes()  { return m_Part.getLocalNodes(); }
  g2l_t getPartLocalCells()  { return m_Part.getLocalCells(); }
  l2l_t getPartN2N()         { return m_Part.getN2N(); }
  l2l_t getPartN2C()         { return m_Part.getN2C(); }
  l2l_t getPartC2C()         { return m_Part.getC2C(); }

public: // methods
  
  Operation();
  virtual ~Operation();
  void del();
    
  void setGrid(vtkUnstructuredGrid *ug) { m_Grid = ug; }
  void setAllCells();
  void setAllVolumeCells();
  void setAllSurfaceCells();
  void setVolume(QString volume_name);
  template <class T> void setCells(const T &cls);
  void setMeshPartition(const MeshPartition &part);

  void setLockGui() { lock_gui = true; }
  OperationThread& getThread() { return thread; }
  void enableAutoSet() { autoset = true; }
  void disableAutoSet() { autoset = false; }
  void setQuickSave(bool b) { m_quicksave = b; }
  void setResetOperationCounter(bool b) { m_resetoperationcounter=b; }
  void setVerboseOn() { m_Verbose = true; }
  void setVerboseOff() { m_Verbose = false; }
  
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

  QString getTypeName() { return m_TypeName; }
  QString getMenuText() { return m_MenuText; }

};

//End of class Operation

Q_DECLARE_INTERFACE(Operation, "eu.engits.enGrid.Operation/1.0")


template <class T>
void Operation::setCells(const T &cls)
{
  m_Part.setGrid(m_Grid);
  m_Part.setCells(cls);
}

#endif
