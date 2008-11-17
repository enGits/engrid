//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#ifndef __vtkEgNormalExtrusion_h
#define __vtkEgNormalExtrusion_h

class vtkEgNormalExtrusion;

#include "vtkEgGridFilter.h"
#include <QVector>

class vtkEgNormalExtrusion : public vtkEgGridFilter
{
  
protected: // attributes
  
  enum mode_t { normal, fixed, planar, cylinder, rotation };
  mode_t mode;
  QVector<double> layer_y;
  vec3_t origin, axis, fixed_normal;
  double min_dist;
  
public: // methods
  
  static vtkEgNormalExtrusion* New();
  void SetLayers(const QVector<double> &y);
  void SetOrigin(vec3_t x) { origin = x; };
  void SetAxis(vec3_t x) { axis = x; };
  void SetNormal(vec3_t x) { fixed_normal = x; };
  void SetMinDist(double d) { min_dist = d; };
  void SetCylindrical() { mode = cylinder; };
  void SetFixed() { mode = fixed; };
  void SetPlanar() { mode = planar; };
  void SetRotation() { mode = rotation; };
  
protected: // methods
  
  vtkEgNormalExtrusion() { mode = normal; };
  ~vtkEgNormalExtrusion() {};
  virtual void ExecuteEg();
  
private: // methods
  
  vtkEgNormalExtrusion (const vtkEgNormalExtrusion&);
  void operator= (const vtkEgNormalExtrusion&);
  
};

#endif
