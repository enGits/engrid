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
#ifndef optimisation_H
#define optimisation_H

class Optimisation;

#include "engrid.h"

class Optimisation
{
  
protected: // attributes
  
  double ***F;
  double Dx;
  double Dy;
  double Dz;
  vec3_t grad_f;
  mat3_t J;
  
protected: // methods
  
  virtual double func(vec3_t x) = 0;
  virtual double func(double x, double y, double z) { return func(vec3_t(x,y,z)); };
  virtual void computeDerivatives(vec3_t x);
  
public: // methods
  
  Optimisation();
   
  vec3_t optimise(vec3_t x);
  void setDeltas(double d) { Dx = d; Dy = d; Dz = d; };
  void setDx(double d) { Dx = d; };
  void setDy(double d) { Dy = d; };
  void setDz(double d) { Dz = d; };
  
};

#endif
