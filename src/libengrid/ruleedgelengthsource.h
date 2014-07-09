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
// 
#ifndef RULEEDGELENGTHSOURCE_H
#define RULEEDGELENGTHSOURCE_H

class RuleEdgeLengthSource;

#include "edgelengthsource.h"
#include "pointfinder.h"

#include <QList>

class RuleEdgeLengthSource : public EdgeLengthSource
{

private: // attributes

  QList<vec3_t> m_Points;
  double        m_GrowthFactor;
  double        m_EdgeLength;
  PointFinder   m_PointFinder;


private: // methods

  void readGrowthFactor();


public:

  RuleEdgeLengthSource(QString rule, vtkUnstructuredGrid *grid);

  virtual double edgeLength(vec3_t x);

};

#endif // RULEEDGELENGTHSOURCE_H
