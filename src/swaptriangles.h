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

#include "surfaceoperation.h"

class SwapTriangles : public SurfaceOperation
{
  
private: // attributes
  
  QVector<bool> m_Marked;
  bool          m_RespectBC;
  bool          m_FeatureSwap;
  int           m_MaxNumLoops;
  
private: // methods
  
  ///returns true if performing a swap on the stencil does not change the orientation of the cells (tetra volume test)
  bool testSwap(stencil_t S);
  
  ///returns true if id_node1 is linked to id_node2
  bool isEdge(vtkIdType id_node1, vtkIdType id_node2);
    
protected: // methods
  
  virtual void operate();

public:

  SwapTriangles();
  void setRespectBC(bool b)   { m_RespectBC   = b; }
  void setFeatureSwap(bool b) { m_FeatureSwap = b; }
  void setMaxNumLoops(int n) { m_MaxNumLoops = n; }
  
};

#endif
