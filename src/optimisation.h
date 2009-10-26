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

  double  m_Err0;
  double  m_ErrS;
  double  m_XS;
  double  m_Exp;
  double  m_MaxErr;
  double  m_TotalError;
  bool    m_Active;
  QString m_Name;
  int     m_NumCalls;

public: // methods

  ErrorFunction();
  void set(QString settings_txt);
  void setName(QString name) { m_Name = name; }
  QString name() { return m_Name; }
  double operator()(double x);
  double maxError() { return m_MaxErr; }
  void reset() { m_MaxErr = 0; m_TotalError = 0; m_NumCalls = 0; }
  bool active() { return m_Active && (m_Err0 > 1e-10); }
  double averageError();
  void activate() { m_Active = true; }
  void deactivate() { m_Active = false; }

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

  void getErrSet(QString group, QString key, double err0, double xs, ErrorFunction &err_func);
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
