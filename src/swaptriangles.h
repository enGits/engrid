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
  
private: // attributes
  
  QVector<bool> m_marked;
  bool m_RespectBC;
  bool m_FeatureSwap;
  
private: // methods
  
  void prepare();
  bool TestSwap(stencil_t S);
  bool isEdge(vtkIdType A, vtkIdType B);
    
protected: // methods
  
  virtual void operate();

public:
  SwapTriangles();
  void setRespectBC(bool b){m_RespectBC=b;};
  void setFeatureSwap(bool b){m_FeatureSwap=b;};
  
};

#endif
