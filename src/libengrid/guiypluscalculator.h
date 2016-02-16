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

#ifndef GUIYPLUSCALCULATOR_H
#define GUIYPLUSCALCULATOR_H

#include "dialogoperation.h"
#include "ui_guiypluscalculator.h"

class GuiYPlusCalculator : public DialogOperation<Ui::GuiYPlusCalculator, Operation>
{

  Q_OBJECT

protected: // attributes

  bool   m_Calculating;
  bool   m_ExternalFlow;
  double m_L;
  double m_U;
  double m_Re;
  double m_Mu;
  double m_Rho;
  double m_Y;
  double m_YPlus;


protected: // methods

  virtual void before();
  virtual void operate() {}

  double cfPlate(double Re_x);
  double cfPipe(double Re);
  void   getValues();
  void   setValues();
  void   calc();


protected slots:

  void changedFlowType(QString);
  void changedL       ();
  void changedU       ();
  void changedMu      ();
  void changedRho     ();
  void changedY       ();
  void changedYPlus   ();
  void changedRe      ();


public:

  GuiYPlusCalculator();

};

#endif // GUIYPLUSCALCULATOR_H
