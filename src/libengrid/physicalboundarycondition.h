// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#ifndef PHYSICALBOUNDARYCONDITION_H
#define PHYSICALBOUNDARYCONDITION_H

#include <QString>
#include <QVector>

#include "engrid.h"

class PhysicalBoundaryCondition
{

private: // attributes

  QString          m_Name;
  QString          m_Type;
  int              m_Index;
  QVector<QString> m_VarNames;
  QVector<double>  m_VarValues;

protected: // methods

public: // methods

  PhysicalBoundaryCondition();

  void setName(QString name) { m_Name = name; }
  void setIndex(int index) { m_Index = index; }
  void setValue(int i, double v) { m_VarValues[i] = v; }
  void setType(QString type);

  QString getName()  { return m_Name; }
  QString getType()  { return m_Type; }
  int     getIndex() { return m_Index; }
  double  getVarValue(int i) { return m_VarValues[i]; }
  QString getVarName(int i)  { return m_VarNames[i]; }
  int     getNumVars()       { return m_VarValues.size(); }

  QString getFoamP();
  QString getFoamU(vec3_t n);
  QString getFoamK();
  QString getFoamEpsilon();
  QString getFoamOmega();
  QString getFoamT();

  QString getFoamType();

};

#endif
