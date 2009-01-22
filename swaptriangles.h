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
#ifndef swaptriangles_H
#define swaptriangles_H

class SwapTriangles;

#include "operation.h"

class SwapTriangles : public Operation
{
  
private: // types
  
  struct stencil_t { 
    vtkIdType id_cell1;
    vtkIdType id_cell2;
    vtkIdType p[4];
    bool valid;
  };
  
private: // attributes
  
  QVector<bool> marked;
  
private: // methods
  
  void prepare();
  stencil_t getStencil(vtkIdType id_cell1, int j1);
    
protected: // methods
  
  virtual void operate();
  
};

#endif
