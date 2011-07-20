// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#ifndef deletecells_H
#define deletecells_H

class DeleteCells;

#include "operation.h"

class DeleteCells : public Operation
{
  
private: // attributes
  
  QVector<vtkIdType> del_cells;
  QVector<vtkIdType> trace_cells;
  
protected: // methods
  
  virtual void operate();
  
public: // methods

  DeleteCells() { EG_TYPENAME; }
  
  template <class T> void setCellsToDelete(const T &cls);
  void setTraceCells(const QVector<vtkIdType> &cells);
  void getTraceCells(QVector<vtkIdType> &cells);
  
};

template <class T>
void DeleteCells::setCellsToDelete(const T &cls)
{
  del_cells.resize(cls.size());
  qCopy(cls.begin(), cls.end(), del_cells.begin());
};


#endif

