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
#ifndef LAPLACESMOOTHER_H
#define LAPLACESMOOTHER_H

#include "operation.h"

class LaplaceSmoother : public Operation {
public:
  LaplaceSmoother();
  virtual void operate();
  bool FlippedCells(vtkIdType id_G, vec3_t P);
    
public:
  void SetBoundaryCodes(QSet<int> a_bcs) { m_bcs=a_bcs; };
  void SetNumberOfIterations(int N){NumberOfIterations=N;};
  
public:
  QSet<int> m_bcs;
  int NumberOfIterations;
  
};

#endif
