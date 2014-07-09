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
#ifndef INSERTPOINTS_H
#define INSERTPOINTS_H

#include "surfaceoperation.h"

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <QSet>
#include <QVector>
#include "egvtkobject.h"
#include "vertexmeshdensity.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

class InsertPoints : public SurfaceOperation
{

private: // data-types

  struct edge_t
  {
    stencil_t S;
    int i_edge;
    double L1, L2, L12;
    double quality() const { return 0.5*min(L1,L2)/L12; }
    bool operator< (const edge_t& E) const { return quality() < E.quality(); }
  };

private: // attributes

  int    m_NumInserted;
  double m_Threshold;

private: // methods
  
  char getNewNodeType(stencil_t S); ///< Returns the type of the node inserted on the edge S.p[1],S.p[3] from stencil_t S

public:

  InsertPoints();
  virtual void operate();
  int insertPoints();
  int getNumInserted() { return m_NumInserted; }
  
};

#endif
