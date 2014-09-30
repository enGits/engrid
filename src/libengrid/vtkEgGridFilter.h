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
#ifndef __vtkEgGridFilter_h
#define __vtkEgGridFilter_h

class vtkEgGridFilter;

#include <vtkUnstructuredGridAlgorithm.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCellLinks.h>
#include <vtkObjectFactory.h>

#include <QVector>

#include "engrid.h"
#include "egvtkobject.h"

class vtkEgGridFilter : public vtkUnstructuredGridAlgorithm, public EgVtkObject
{
  
protected: // attributes
  
  QSet<int>            m_BoundaryCodes;
  vtkUnstructuredGrid *m_Input;
  vtkUnstructuredGrid *m_Output;
  vtkUnstructuredGrid *m_LastOutput;
  unsigned long        m_LastRun;
  
protected: // methods
  
  vtkEgGridFilter();
  ~vtkEgGridFilter();
  
  virtual int FillInputPortInformation(int, vtkInformation* info);
  virtual int FillOutputPortInformation(int, vtkInformation* info);
  virtual int RequestData(vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector*  outputVector);
  
  virtual void ExecuteEg() = 0;
  
  /**
   * Extract boundary elements and nodes with a given set of boundary conditions
   * @param cells On return, this will contain the indexes of the extracted boundary elements.
   * @param modes On return, this will contain the indexes of the extracted boundary nodes.
   * @param grid  The grid to operate on.
   * @param bc_c  A Qt container with the boundary conditions to extract.
   */
  void ExtractBoundary(QVector<vtkIdType>  &cells, QVector<vtkIdType> &nodes, QSet<int> &bc_c, vtkUnstructuredGrid *grid);
  
public: // methods
  
  template <typename C>
  void SetBoundaryCodes(const C& bc);
  virtual void Update() { vtkAlgorithm::Update(); ExecuteEg(); }

};

template <typename C>
void vtkEgGridFilter::SetBoundaryCodes(const C &bc_c)
{
  bool update = false;

  QSet<int> bc_set;
  foreach (int bc , bc_c) {
    bc_set.insert(bc);
  }

  if (m_BoundaryCodes.size() != bc_set.size()) {
    update = true;
  } else {
    QSet<int> bc_inters = m_BoundaryCodes.intersect(bc_set);
    if (bc_inters.size() != bc_set.size()) {
      update = true;
    }
  }
  if (update) {
    m_BoundaryCodes = bc_set;
    Modified();
  }
}


#endif
