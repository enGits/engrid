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
#ifndef FIXCADGEOMETRY_H
#define FIXCADGEOMETRY_H

#include "surfacealgorithm.h"

class FixCadGeometry: public SurfaceAlgorithm
{

private: // attributes

  int    m_NumNonManifold;
  double m_OriginalFeatureAngle;
  
protected: // methods
  
  void customUpdateNodeInfo();
  void callMesher();
  void setDesiredLength(double L = 1e99);
  void copyFaces(const QVector<bool> &copy_face);
  void fixNonManifold1();
  void fixNonManifold2();
  void markNonManifold();

  virtual void operate();
  
public: // methods
  
  FixCadGeometry();
  
};

#endif
