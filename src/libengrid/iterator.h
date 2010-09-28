//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#ifndef iterator_H
#define iterator_H

class Iterator;

#include "operation.h"

class Iterator : public Operation
{
  
private: // attributes
  
protected: // attributes
  
  struct pair_t 
  {
    vtkIdType item1;
    vtkIdType item2;
    bool      terminate;
  };
  
  QVector<pair_t> pair;
  
  /** array used to mark item for the next iteration loop. */
  QVector<bool> mark1;
  
  /** array used to marked "finished" items */
  QVector<bool> mark2;
  
  QVector<vtkIdType>  item;
  
  bool volume_iteration;
  bool custom_iteration;
  
  QVector<vtkIdType>     cells;
  QVector<vtkIdType>     nodes;
  QVector<int>           _cells;
  QVector<int>           _nodes;
  QVector<QSet<int> >    n2c;
  QVector<QSet<int> >    n2n;
  QVector<QVector<int> > c2c;
  
  
protected: // methods
  
  virtual void pass1();
  virtual void pass2() {}
  void getCells();
  
public: // methods
  
  Iterator();
  
  void setVolumeIteration()  { volume_iteration = true; }
  void setSurfaceIteration() { volume_iteration = false; }
  
};

#endif

