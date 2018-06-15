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
#include "cgaltricadinterface.h"

class FixCadGeometry: public SurfaceAlgorithm
{

private: // attributes

  int    m_NumNonManifold;
  double m_OriginalFeatureAngle;
  double m_SnapTolerance;

  CgalTriCadInterface* m_Cad;


protected: // data types

  struct cut_t
  {
    vec3_t x;
    double w;
    double R;
    double L;
    bool   edge_cut;
    bool   node1_surf;
    bool   node2_surf;
  };


  
protected: // methods
  
  void customUpdateNodeInfo();
  void callMesher();
  void setDesiredLength(double L = 1e99);
  void copyFaces(const QVector<bool> &copy_face);
  void fixNonManifold1();
  void fixNonManifold2();
  void markNonManifold();

  void  computeCharLength();
  void  createBox();
  cut_t snapCut(vtkIdType id_node1, vtkIdType id_node2);
  void  marchOutside();
  void  cut();
  void  refine();

  virtual void operate();
  
public: // methods
  
  FixCadGeometry();

  void setSnapTolerance(double tol) { m_SnapTolerance = tol; }
  
};

#endif
