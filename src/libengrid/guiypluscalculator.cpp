// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2016 enGits GmbH                                      +
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

#include "guiypluscalculator.h"

GuiYPlusCalculator::GuiYPlusCalculator()
{
  m_Ui.m_LineEditL->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditRe->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditMu->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditU->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditRho->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditY->setValidator(new QDoubleValidator());
  m_Ui.m_LineEditYPlus->setValidator(new QDoubleValidator());
  connect(m_Ui.m_LineEditL,        SIGNAL(editingFinished()),            this, SLOT(changedL()));
  connect(m_Ui.m_LineEditMu,       SIGNAL(editingFinished()),            this, SLOT(changedMu()));
  connect(m_Ui.m_LineEditRe,       SIGNAL(editingFinished()),            this, SLOT(changedRe()));
  connect(m_Ui.m_LineEditRho,      SIGNAL(editingFinished()),            this, SLOT(changedRho()));
  connect(m_Ui.m_LineEditU,        SIGNAL(editingFinished()),            this, SLOT(changedU()));
  connect(m_Ui.m_LineEditY,        SIGNAL(editingFinished()),            this, SLOT(changedY()));
  connect(m_Ui.m_LineEditYPlus,    SIGNAL(editingFinished()),            this, SLOT(changedYPlus()));
  connect(m_Ui.m_ComboBoxFlowType, SIGNAL(currentIndexChanged(QString)), this, SLOT(changedFlowType(QString)));
}

double GuiYPlusCalculator::cfPlate(double Re_x)
{
  return 0.455/pow(log(0.06*Re_x), 2.0);
}

double GuiYPlusCalculator::cfPipe(double Re)
{
  double lam = 0.01;
  for (int iter = 0; iter < 100; ++iter) {
    lam = 1.0/pow((2*log10(Re*sqrt(lam)) - 0.8), 2.0);
  }
  return 0.25*lam;
}

void GuiYPlusCalculator::before()
{
  m_Calculating = false;
  m_ExternalFlow = false;
  getValues();
  calc();
}

void GuiYPlusCalculator::getValues()
{
  m_Calculating = true;
  m_L     = m_Ui.m_LineEditL->text().toDouble();
  m_U     = m_Ui.m_LineEditU->text().toDouble();
  m_Re    = m_Ui.m_LineEditRe->text().toDouble();
  m_Mu    = m_Ui.m_LineEditMu->text().toDouble();
  m_Rho   = m_Ui.m_LineEditRho->text().toDouble();
  m_Y     = m_Ui.m_LineEditY->text().toDouble();
  m_YPlus = m_Ui.m_LineEditYPlus->text().toDouble();
  m_ExternalFlow = m_Ui.m_ComboBoxFlowType->currentText() == "external";
}

void GuiYPlusCalculator::setValues()
{
  QString num;
  num.setNum(m_L);     m_Ui.m_LineEditL->setText(num);
  num.setNum(m_U);     m_Ui.m_LineEditU->setText(num);
  num.setNum(m_Re);    m_Ui.m_LineEditRe->setText(num);
  num.setNum(m_Mu);    m_Ui.m_LineEditMu->setText(num);
  num.setNum(m_Rho);   m_Ui.m_LineEditRho->setText(num);
  num.setNum(m_Y);     m_Ui.m_LineEditY->setText(num);
  num.setNum(m_YPlus); m_Ui.m_LineEditYPlus->setText(num);
  m_Calculating = false;
}

void GuiYPlusCalculator::calc()
{
  double nu = m_Mu/m_Rho;
  m_Re = m_U*m_L/nu;
  double cf = cfPlate(m_Re);
  if (!m_ExternalFlow) {
    cf = cfPipe(m_Re);
  }
  m_YPlus = m_U*m_Y/nu*sqrt(0.5*cf);
  setValues();
}

void GuiYPlusCalculator::changedFlowType(QString)
{
  if (m_Calculating) {
    return;
  }
  getValues();
  calc();
}

void GuiYPlusCalculator::changedL()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  calc();
}

void GuiYPlusCalculator::changedMu()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  calc();
}

void GuiYPlusCalculator::changedRe()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  m_Re = m_Ui.m_LineEditRe->text().toDouble();
  m_U = m_Re*m_Mu/(m_Rho*m_L);
  calc();
}

void GuiYPlusCalculator::changedRho()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  m_Rho = m_Ui.m_LineEditRho->text().toDouble();
  calc();
}

void GuiYPlusCalculator::changedU()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  calc();
}

void GuiYPlusCalculator::changedY()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  calc();
}

void GuiYPlusCalculator::changedYPlus()
{
  if (m_Calculating) {
    return;
  }
  getValues();
  double nu = m_Mu/m_Rho;
  m_Re = m_U*m_L/nu;
  double cf = cfPlate(m_Re);
  if (!m_ExternalFlow) {
    cf = cfPipe(m_Re);
  }
  m_YPlus = m_Ui.m_LineEditYPlus->text().toDouble();
  m_Y = sqrt(2.0/cf)*m_YPlus*nu/m_U;
  calc();
}
