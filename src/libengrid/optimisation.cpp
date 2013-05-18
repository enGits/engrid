// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
  m_Err0 = 1.0;
  m_ErrS = 0.5;
  m_XS   = 0.5;
  m_Exp  = 1.0;

  m_TotalError = 0.0;
  m_Active     = true;
}

void ErrorFunction::set(QString settings_txt)
{
  QStringList items = settings_txt.split(';');
  if (items.size() != 2) {
    EG_ERR_RETURN("syntax error for error weighting");
  }
  m_Err0 = items[0].trimmed().toDouble();
  m_ErrS = 0.5*m_Err0;
  m_XS   = items[1].trimmed().toDouble();
  m_Exp = 1.0;
  if (m_Err0 > 1e-3) {
    m_Exp = log(m_ErrS/m_Err0)/log(1.0 - m_XS);
  }
}

double ErrorFunction::operator()(double x)
{
  double e = 0;
  if (m_Active) {
    double x0 = fabs(1-x);
    e = m_Err0*pow(x0, m_Exp);
    m_MaxErr = max(m_MaxErr, x0);
    m_TotalError += e;
    m_AverageError += x0;
    ++m_NumCalls;
    ++m_NumTotalCalls;
  }
  return e;
}

double ErrorFunction::averageError()
{
  if (m_NumTotalCalls == 0) {
    return 0;
  }
  return m_AverageError/m_NumTotalCalls;
}

double ErrorFunction::totalError()
{
  if (m_NumCalls == 0) {
    return 0;
  }
  return m_TotalError/m_NumCalls;
}

void ErrorFunction::reset(bool reset_average)
{
  m_MaxErr = 0;
  m_TotalError = 0;
  m_NumCalls = 0;
  if (reset_average) {
    m_AverageError = 0;
    m_NumTotalCalls = 0;
  }
}

Optimisation::Optimisation()
{
  setDeltas(1e-6);
  F = new double**[3];
  for (int i = 0; i < 3; ++i) {
    F[i] = new double*[3];
    for (int j = 0; j < 3; ++j) {
      F[i][j] = new double[3];
    }
  }
  m_ErrorFunctions.clear();
}

void Optimisation::getErrSet(QString group, QString key, double err0, double xs, ErrorFunction* err_func)
{
  QString err0_txt, xs_txt;
  err0_txt.setNum(err0);
  xs_txt.setNum(xs);
  QString value = err0_txt + "; " + xs_txt;
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
  err_func->set(variable);
  err_func->setName(key.replace(" ", "_"));
  m_ErrorFunctions.append(err_func);
}

double Optimisation::angleX(const vec3_t &v1, const vec3_t &v2)
{
  double va1  = v1.abs();
  double va2  = v2.abs();
  double scal = max(0.0, v1*v2);
  double alpha = acos(scal/(va1*va2));
  return fabs(1 - 2*alpha/M_PI);
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
}

vec3_t Optimisation::optimise(vec3_t x) 
{ 
  computeDerivatives(x);
  mat3_t M = J.inverse();
  return (-1)*(M*grad_f);
}


void Optimisation::resetErrorFunctions(bool reset_average)
{
  foreach (ErrorFunction *err_func, m_ErrorFunctions) {
    err_func->reset(reset_average);
  }
}

double Optimisation::totalError()
{
  double e = 0;
  foreach (ErrorFunction *err_func, m_ErrorFunctions) {
    e += err_func->totalError();
  }
  return e;
}

void Optimisation::printErrors()
{
  foreach (ErrorFunction *err_func, m_ErrorFunctions) {
    if (err_func->active()) {
      cout << qPrintable(err_func->name()) << ":\n  average=" << err_func->averageError() << "\n  maximum=" << err_func->maxError() << endl;
    }
  }
}
