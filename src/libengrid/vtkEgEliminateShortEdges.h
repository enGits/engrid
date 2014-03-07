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
#ifndef __vtkEgEliminateShortEdges_h
#define __vtkEgEliminateShortEdges_h

class vtkEgEliminateShortEdges;

#include "vtkEgGridFilter.h"

#include <QVector>

/**
 * A VTK filter to remove elements with bad aspect ratios.
 * Currently this will only work with triangular surface grids.
 * There is also no error checking done -- yet -- which means
 * that the resulting grid has to be carefully looked at.
 * This filter has, however, proven to be useful in certain
 * cases.
 */
class vtkEgEliminateShortEdges : public vtkEgGridFilter
{
  
private: // attributes
  
  QVector<vtkIdType> new_node_idx;
  QVector<vtkIdType> node_mapping;
  QVector<bool> delete_cell;
  double max_ratio;
  double max_length;
  int N_eliminated;
  vtkIdType N_new_points;
  vtkIdType N_new_cells;
  
public: // static methods
  
  static vtkEgEliminateShortEdges* New();
  
public: // methods
  
  void SetMaxRatio(double mr)  { max_ratio  = mr ; Modified(); };
  void SetMaxLength(double ml) { max_length = ml ; Modified(); };
  
  int getNumEliminated() { return N_eliminated; };
  
protected: // methods
  
  vtkEgEliminateShortEdges();
  ~vtkEgEliminateShortEdges() {};
  virtual void ExecuteEg();
  
private: // methods
  
  vtkEgEliminateShortEdges (const vtkEgEliminateShortEdges&);
  void operator= (const vtkEgEliminateShortEdges&);
  
  void CheckEdges();
  void CheckCells();
  void CopyPoints();
  void CopyCells();
  
};

#endif
