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
#ifndef operation_H
#define operation_H

class Operation;
class GuiMainWindow;

#include "egvtkobject.h"

#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkSmartPointer.h>

#include <QThread>
#include <QMutex>

#include <typeinfo>

class OperationThread : public QThread
{
  
private:
  
  Operation *op;
  
protected:
  
  virtual void run();
  
public:
  
  void setOperation(Operation *an_op) { op = an_op; };
  
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
  
  QVector<vtkIdType> nodes_map;
  QVector<vtkIdType> cells_map;
  bool               gui;
  bool               autoset;
  Error             *err;
  
private: // methods
  
  void initMapping();
  
protected: // attributes
  
  vtkUnstructuredGrid   *grid;
  QVector<vtkIdType>     cells;
  QVector<int>           _cells;
  QVector<vtkIdType>     nodes;
  QVector<int>           _nodes;
  QVector<QSet<int> >    n2c;
  QVector<QSet<int> >    n2n;
  QVector<QVector<int> > c2c;
  QVector<bool>          node_fixed;
  QVector<bool>          cell_fixed;
  
protected: // methods
  
  void checkGrid();
  void updateActors();
  GuiMainWindow* mainWindow();
  virtual void operate() = 0;

public: // methods
  
  Operation();
  virtual ~Operation();
  void del() { garbage_operations.insert(this); };
    
  void setGrid(vtkUnstructuredGrid *ug) { grid = ug; };
  void setAllCells();
  void setAllVolumeCells();
  void setAllSurfaceCells();
  vtkIdType getNewNode(vtkIdType id_old_node) { return nodes_map[_nodes[id_old_node]] ; };
  vtkIdType getNewCell(vtkIdType id_old_cell) { return cells_map[_cells[id_old_cell]] ; };
  void setNewNode(vtkIdType id_old_node, vtkIdType id_new_node) { nodes_map[_nodes[id_old_node]] = id_new_node; };
  void setNewCell(vtkIdType id_old_cell, vtkIdType id_new_cell) { cells_map[_cells[id_old_cell]] = id_new_cell; };
  void setGui() { gui = true; };
  OperationThread& getThread() { return thread; };
  void enableAutoSet() { autoset = true; };
  void disableAutoSet() { autoset = false; };
  
  virtual void operator()();

  template <class T> void setCells(const T &cls);
  template <class T> void setNodes(const T &nds);

  static void collectGarbage();
  
};

template <class T>
void Operation::setCells(const T &cls)
{
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  getNodesFromCells(cells, nodes, grid);
  createCellMapping(cells, _cells, grid);
  createNodeMapping(nodes, _nodes, grid);
  createNodeToCell(cells, nodes, _nodes, n2c, grid);
  createNodeToNode(cells, nodes, _nodes, n2n, grid);
  createCellToCell(cells, c2c, grid);
  node_fixed.fill(nodes.size(), false);
  cell_fixed.fill(cells.size(), false);
  initMapping();
};

template <class T>
void Operation::setNodes(const T &nds)
{
  nodes.resize(nds.size());
  qCopy(nds.begin(), nds.end(), nodes.begin());
  createNodeMapping(nodes, _nodes, grid);
  QSet<vtkIdType> cls;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType *pts, N_pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
      if (_nodes[pts[i_pts]] >= 0) {
        cls.insert(id_cell);
        break;
      };
    };
  };
  cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), cells.begin());
  createCellMapping(cells, _cells, grid);
  createNodeToCell(cells, nodes, _nodes, n2c, grid);
  createNodeToNode(cells, nodes, _nodes, n2n, grid);
  createCellToCell(cells, c2c, grid);
  node_fixed.fill(nodes.size(), false);
  cell_fixed.fill(cells.size(), false);
  initMapping();
};


#endif
