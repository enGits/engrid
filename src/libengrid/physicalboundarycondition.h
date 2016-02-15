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
  QVector<QString> m_VarValues;
  QVector<QString> m_VarTypes;

protected: // methods

  void addBoolVar  (QString name, bool v);
  void addRealVar  (QString name, double v);
  void addStringVar(QString name, QString v);
  void addIntVar   (QString name, int v);
  void addVecVar   (QString name, vec3_t v);


public: // methods

  PhysicalBoundaryCondition();

  void setName(QString name)      { m_Name = name; }
  void setIndex(int index)        { m_Index = index; }

  void checkVarType(int i, QString type);
  void setValue(int i, double v);
  void setValue(int i, QString v);
  void setValue(int i, int v);
  void setValue(int i, bool v);
  void setValue(int i, vec3_t v);
  void setValueFromString(int i, QString v);
  void setType(QString type);

  QString getName()  { return m_Name; }
  QString getType()  { return m_Type; }
  int     getIndex() { return m_Index; }
  QString getVarValueAsString(int i) { return m_VarValues[i]; }
  int     getVarValueAsInt(int i)    { return m_VarValues[i].toInt(); }
  double  getVarValueAsDouble(int i) { return m_VarValues[i].toDouble(); }
  bool    getVarValueAsBool(int i);
  vec3_t  getVarValueAsVec3(int i);
  QString getVarType(int i)          { return m_VarTypes[i]; }
  QString getVarName(int i)          { return m_VarNames[i]; }
  int     getNumVars()               { return m_VarValues.size(); }

  QString getFoamP(QString version);
  QString getFoamU(QString version, vec3_t n);
  QString getFoamK(QString version);
  QString getFoamEpsilon(QString version);
  QString getFoamOmega(QString version);
  QString getFoamT(QString version);
  QString getFoamNut(QString version);
  QString getFoamType();

  QString xmlText();

};

#endif
