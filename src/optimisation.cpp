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
#include "optimisation.h"
#include "guimainwindow.h"

ErrorFunction::ErrorFunction()
{
  m_Weighting1 = 1.0;
  m_Weighting2 = 1.0;
  m_XSwitch = 0.5;
  m_Exponent = 1.0;
  m_TotalError = 0.0;
  m_Active = true;
}

void ErrorFunction::set(QString settings_txt)
{
  QStringList items = settings_txt.split(',');
  if (items.size() != 4) {
    EG_ERR_RETURN("syntax error for error weighting");
  }
  m_Weighting1 = items[0].trimmed().toDouble();
  m_Weighting2 = items[1].trimmed().toDouble();
  m_Exponent   = items[2].trimmed().toDouble();
  m_XSwitch    = items[3].trimmed().toDouble();
}

double ErrorFunction::operator ()(double x)
{
  double e1 = max(0.0, -m_Weighting1*(x - m_XSwitch));
  double e2 = fabs(1 - x);
  double e  = max(e1, m_Weighting2*pow(e2, m_Exponent));
  m_MaxErr = max(m_MaxErr, e2);
  m_TotalError += e2;
  ++m_NumCalls;
  return e;
}

double ErrorFunction::averageError()
{
  if (m_NumCalls == 0) {
    return 0;
  }
  return m_TotalError/m_NumCalls;
}


Optimisation::Optimisation()
{
  setDeltas(1e-6);
  F = new double**[3];
  for (int i = 0; i < 3; ++i) {
    F[i] = new double*[3];
    for (int j = 0; j < 3; ++j) {
      F[i][j] = new double[3];
    };
  };
};

void Optimisation::getErrSet(QString group, QString key, double w1, double w2, double e, double s, ErrorFunction &err_func)
{
  QString w1_txt, w2_txt, s_txt, e_txt;
  w1_txt.setNum(w1);
  w2_txt.setNum(w2);
  s_txt.setNum(s);
  e_txt.setNum(e);
  QString value = w1_txt + ", " + w2_txt + ", " + e_txt + ", " + s_txt;
  QString variable;
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "string/" + key;
  if (group != QObject::tr("General")) {
    qset->beginGroup(group);
  }
  if (!qset->contains(typed_key)) {
    qset->setValue(typed_key, value);
  }
  variable = (qset->value(typed_key,variable)).toString();
  if (group != QObject::tr("General")) {
    qset->endGroup();
  }
  err_func.set(variable);
  err_func.setName(key.replace(" ", "_"));
}

double Optimisation::angleX(const vec3_t &v1, const vec3_t &v2)
{
  double scal = v1*v2;
  double alpha = acos(scal);
  return fabs(1 - alpha/M_PI);
}

void Optimisation::computeDerivatives(vec3_t x)
{
  F[0][0][0] = func(x[0]-Dx, x[1]-Dy, x[2]-Dz);
  F[1][0][0] = func(x[0],    x[1]-Dy, x[2]-Dz);
  F[2][0][0] = func(x[0]+Dx, x[1]-Dy, x[2]-Dz);
  F[0][1][0] = func(x[0]-Dx, x[1],    x[2]-Dz);
  F[1][1][0] = func(x[0],    x[1],    x[2]-Dz);
  F[2][1][0] = func(x[0]+Dx, x[1],    x[2]-Dz);
  F[0][2][0] = func(x[0]-Dx, x[1]+Dy, x[2]-Dz);
  F[1][2][0] = func(x[0],    x[1]+Dy, x[2]-Dz);
  F[2][2][0] = func(x[0]+Dx, x[1]+Dy, x[2]-Dz);
  F[0][0][1] = func(x[0]-Dx, x[1]-Dy, x[2]);
  F[1][0][1] = func(x[0],    x[1]-Dy, x[2]);
  F[2][0][1] = func(x[0]+Dx, x[1]-Dy, x[2]);
  F[0][1][1] = func(x[0]-Dx, x[1],    x[2]);
  F[1][1][1] = func(x[0],    x[1],    x[2]);
  F[2][1][1] = func(x[0]+Dx, x[1],    x[2]);
  F[0][2][1] = func(x[0]-Dx, x[1]+Dy, x[2]);
  F[1][2][1] = func(x[0],    x[1]+Dy, x[2]);
  F[2][2][1] = func(x[0]+Dx, x[1]+Dy, x[2]);
  F[0][0][2] = func(x[0]-Dx, x[1]-Dy, x[2]+Dz);
  F[1][0][2] = func(x[0],    x[1]-Dy, x[2]+Dz);
  F[2][0][2] = func(x[0]+Dx, x[1]-Dy, x[2]+Dz);
  F[0][1][2] = func(x[0]-Dx, x[1],    x[2]+Dz);
  F[1][1][2] = func(x[0],    x[1],    x[2]+Dz);
  F[2][1][2] = func(x[0]+Dx, x[1],    x[2]+Dz);
  F[0][2][2] = func(x[0]-Dx, x[1]+Dy, x[2]+Dz);
  F[1][2][2] = func(x[0],    x[1]+Dy, x[2]+Dz);
  F[2][2][2] = func(x[0]+Dx, x[1]+Dy, x[2]+Dz);
  
  grad_f[0] = (F[2][1][1]-F[0][1][1])/(2*Dx);
  grad_f[1] = (F[1][2][1]-F[1][0][1])/(2*Dy);
  grad_f[2] = (F[1][1][2]-F[1][1][0])/(2*Dz);
  
  J[0][0] = (F[0][1][1]-2*F[1][1][1]+F[2][1][1])/(Dx*Dx);
  J[1][1] = (F[1][0][1]-2*F[1][1][1]+F[1][2][1])/(Dy*Dy);
  J[2][2] = (F[1][1][0]-2*F[1][1][1]+F[1][1][2])/(Dz*Dz);
  
  J[0][1] = ((F[2][2][1]-F[0][2][1]) - (F[2][0][1]-F[0][0][1]))/(4*Dx*Dy);
  J[0][2] = ((F[2][1][2]-F[0][1][2]) - (F[2][1][0]-F[0][1][0]))/(4*Dx*Dz);
  J[1][2] = ((F[1][2][2]-F[1][2][2]) - (F[1][0][0]-F[1][0][0]))/(4*Dy*Dz);
  J[1][0] = J[0][1];
  J[2][0] = J[0][2];
  J[2][1] = J[1][2];
};

vec3_t Optimisation::optimise(vec3_t x) 
{ 
  computeDerivatives(x);
  mat3_t M = J.inverse();
  return (-1)*(M*grad_f);
};
