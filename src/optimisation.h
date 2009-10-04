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
#ifndef optimisation_H
#define optimisation_H

class Optimisation;
class ErrorFunction;

#include "engrid.h"

class ErrorFunction
{

private: // attributes

  double m_Weighting1;
  double m_Weighting2;
  double m_XSwitch;
  double m_Exponent;
  double m_MaxErr;
  double m_LastError;

public: // methods

  ErrorFunction();
  void set(QString settings_txt);
  double operator()(double x);
  double maxError() { return m_MaxErr; }
  void reset() { m_MaxErr = 0; }
  bool active() { return (m_Weighting1 > 1e-10) || (m_Weighting2 > 1e-10); }
  double lastError() { return m_LastError; }

};


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

  void getErrSet(QString group, QString key, double w1, double w2, double e, double s, ErrorFunction &err_func);
  double angleX(const vec3_t &v1, const vec3_t &v2);
  
public: // methods
  
  Optimisation();
   
  vec3_t optimise(vec3_t x);
  void setDeltas(double d) { Dx = d; Dy = d; Dz = d; };
  void setDx(double d) { Dx = d; };
  void setDy(double d) { Dy = d; };
  void setDz(double d) { Dz = d; };
  
};

#endif
