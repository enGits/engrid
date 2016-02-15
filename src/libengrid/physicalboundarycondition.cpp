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
#include "physicalboundarycondition.h"

#include <QStringList>
#include <QtDebug>

PhysicalBoundaryCondition::PhysicalBoundaryCondition()
{
  m_Name = "unknown";
  m_Index = -1;
  m_Type = "unknown";
}

void PhysicalBoundaryCondition::addBoolVar(QString name, bool v)
{
  m_VarTypes.push_back("bool");
  m_VarNames.push_back(name);
  if (v) {
    m_VarValues.push_back("yes");
  } else {
    m_VarValues.push_back("no");
  }
}

void PhysicalBoundaryCondition::addIntVar(QString name, int v)
{
  m_VarTypes.push_back("int");
  m_VarNames.push_back(name);
  QString txt;
  txt.setNum(v);
  m_VarValues.push_back(txt);
}

void PhysicalBoundaryCondition::addRealVar(QString name, double v)
{
  m_VarTypes.push_back("real");
  m_VarNames.push_back(name);
  QString txt;
  txt.setNum(v);
  m_VarValues.push_back(txt);
}

void PhysicalBoundaryCondition::addStringVar(QString name, QString v)
{
  m_VarTypes.push_back("string");
  m_VarNames.push_back(name);
  m_VarValues.push_back(v);
}

void PhysicalBoundaryCondition::addVecVar(QString name, vec3_t v)
{
  m_VarTypes.push_back("vector");
  m_VarNames.push_back(name);
  QString txt, num;
  txt = "(";
  num.setNum(v[0]);
  txt += num;
  txt += ", ";
  num.setNum(v[1]);
  txt += num;
  txt += ", ";
  num.setNum(v[2]);
  txt += num;
  txt += ")";
  m_VarValues.push_back(txt);
}

void PhysicalBoundaryCondition::checkVarType(int i, QString type)
{
  if (i >= m_VarTypes.size())  EG_BUG;
  if (i >= m_VarNames.size())  EG_BUG;
  if (i >= m_VarValues.size()) EG_BUG;
  if (m_VarTypes[i] != type)   EG_BUG;
}

void PhysicalBoundaryCondition::setValue(int i, bool v)
{
  checkVarType(i, "bool");
  if (v) {
    m_VarValues[i] = "yes";
  } else {
    m_VarValues[i] = "no";
  }
}

void PhysicalBoundaryCondition::setValue(int i, double v)
{
  checkVarType(i, "real");
  QString num;
  num.setNum(v);
  m_VarValues[i] = num;
}

void PhysicalBoundaryCondition::setValue(int i, int v)
{
  checkVarType(i, "int");
  QString num;
  num.setNum(v);
  m_VarValues[i] = num;
}

void PhysicalBoundaryCondition::setValue(int i, QString v)
{
  checkVarType(i, "string");
  m_VarValues[i] = v;
}

void PhysicalBoundaryCondition::setValue(int i, vec3_t v)
{
  checkVarType(i, "vector");
  QString num, txt = "(";
  num.setNum(v[0]);
  txt += num;
  txt += ", ";
  num.setNum(v[1]);
  txt += num;
  txt += ", ";
  num.setNum(v[2]);
  txt += num;
  txt += ")";
  m_VarValues[i] = txt;
}

void PhysicalBoundaryCondition::setValueFromString(int i, QString v)
{
  if (i >= m_VarTypes.size())  EG_BUG;
  if (i >= m_VarNames.size())  EG_BUG;
  if (i >= m_VarValues.size()) EG_BUG;
  m_VarValues[i] = v;
  if (getVarType(i) == "real") {
    setValue(i, getVarValueAsDouble(i));
  }
  if (getVarType(i) == "int") {
    setValue(i, getVarValueAsInt(i));
  }
  if (getVarType(i) == "bool") {
    setValue(i, getVarValueAsBool(i));
  }
  if (getVarType(i) == "vector") {
    setValue(i, getVarValueAsVec3(i));
  }
}

QString PhysicalBoundaryCondition::xmlText()
{
  QString txt;
  txt.setNum(m_Index);
  txt += " " + m_Name + " " + m_Type + ";";
  for (int i = 0; i < m_VarTypes.size(); ++i) {
    txt += " " + m_VarTypes[i] + " " + m_VarNames[i] + " = " + m_VarValues[i] + ";";
  }
  return txt;
}

bool PhysicalBoundaryCondition::getVarValueAsBool(int i)
{
  QString value = m_VarValues[i].toLower();
  if (value == "no")    return false;
  if (value == "off")   return false;
  if (value == "false") return false;
  if (value == "yes")   return true;
  if (value == "on")    return true;
  if (value == "true")  return true;
  EG_ERR_RETURN("boundary condition \"" + m_Name + "\": cannot convert \"" + m_VarValues[i] + "\" to boolen (true/false)!");
  return false;
}

vec3_t PhysicalBoundaryCondition::getVarValueAsVec3(int i)
{
  QString txt = m_VarValues[i];
  txt = txt.replace("(", "");
  txt = txt.replace(")", "");
  txt = txt.trimmed();
  QStringList parts = txt.split(",");
  vec3_t v(0,0,0);
  if (parts.size() == 3) {
    for (int j = 0; j < 3; ++j) {
      v[j] = parts[j].toDouble();
    }
  }
  return v;
}

void PhysicalBoundaryCondition::setType(QString type)
{
  m_Type = type;
  m_VarNames.clear();
  m_VarValues.clear();
  if (m_Type == "turbulent-duct-inlet") {
    addRealVar("velocity", 1);
    addRealVar("temperature", 300);
  }
  if (m_Type == "laminar-duct-inlet") {
    addRealVar("velocity", 1);
    addRealVar("temperature", 300);
  }
  if (m_Type == "inflow") {
    addRealVar("velocity", 1);
    addRealVar("turbulent-intensity", 0.04);
    addRealVar("turbulent-length-scale", 1);
    addRealVar("temperature", 300);
  }
  if (m_Type == "outflow") {
    addBoolVar("supersonic", false);
    addRealVar("pressure", 0);
  }
  if (m_Type == "cyclic") {
    addStringVar("counterpart", "");
  }
  if (m_Type == "symmetry") {
  }
  if (m_Type == "DrNUM-turbulent-wall") {
    addRealVar("yplus-cut-off", 300);
  }
  if (m_Type == "wall") {
  }
  if (m_Type == "inviscid wall") {
  }
}

QString PhysicalBoundaryCondition::getFoamEpsilon(QString version)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    if (version >= "2.1") {
      s << "        type    omegaWallFunction;\n";
      s << "        Cmu     0.09;\n";
      s << "        kappa   0.41;\n";
      s << "        E       9.8;\n";
      s << "        beta1   0.075;\n";
      s << "        value   uniform 0;\n";
    } else {
      s << "        type zeroGradient;\n";
    }
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    double k       = 1.5*sqr(getVarValueAsDouble(0)*getVarValueAsDouble(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValueAsDouble(2);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << epsilon << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamK(QString version)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    if (version >= "2.1") {
      s << "        type    kqRWallFunction;\n";
      s << "        value   uniform 0;\n";
    } else {
      s << "        type zeroGradient;\n";
    }
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    double k       = 1.5*sqr(getVarValueAsDouble(0)*getVarValueAsDouble(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValueAsDouble(2);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << k << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamOmega(QString version)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    if (version >= "2.1") {
      s << "        type    omegaWallFunction;\n";
      s << "        Cmu     0.09;\n";
      s << "        kappa   0.41;\n";
      s << "        E       9.8;\n";
      s << "        beta1   0.075;\n";
      s << "        value   uniform 0;\n";
    } else {
      s << "        type zeroGradient;\n";
    }
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    double k       = 1.5*sqr(getVarValueAsDouble(0)*getVarValueAsDouble(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValueAsDouble(2);
    double omega   = epsilon/(0.09*k);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << omega << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamNut(QString version)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    if (version >= "2.1") {
      s << "        type    nutkWallFunction;\n";
      s << "        Cmu     0.09;\n";
      s << "        kappa   0.41;\n";
      s << "        E       9.8;\n";
      s << "        value   uniform 0;\n";
    } else {
      EG_BUG;
    }
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    double k       = 1.5*sqr(getVarValueAsDouble(0)*getVarValueAsDouble(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValueAsDouble(2);
    double omega   = epsilon/(0.09*k);
    s << "        type  calculated;\n";
    s << "        value uniform 0;\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamP(QString)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "outlet") {
    s << "        type  fixedValue;\n";
    s << "        value uniform " << getVarValueAsDouble(0) << ";\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamU(QString, vec3_t n)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    s << "        type  fixedValue;\n";
    s << "        value uniform (0 0 0);\n";
  }
  if (m_Type == "slip") {
    s << "        type slip;\n";
  }
  if (m_Type == "inlet") {
    s << "        type  fixedValue;\n";
    s << "        value uniform (" << getVarValueAsDouble(0)*n[0] << " " << getVarValueAsDouble(0)*n[1] << " " << getVarValueAsDouble(0)*n[2] << ");\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamT(QString version)
{
  QString str;
  QTextStream s(&str, QIODevice::WriteOnly);
  if (m_Type == "symmetry") {
    s << "        type symmetryPlane;\n";
  }
  if (m_Type == "wall") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "slip") {
    s << "        type zeroGradient;\n";
  }
  if (m_Type == "inlet") {
    s << "        type  fixedValue;\n";
    s << "        value uniform " << getVarValueAsDouble(3) << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamType()
{
  if (m_Type == "symmetry") {
    return ("symmetryPlane");
  }
  if (m_Type == "wall") {
    return ("wall");
  }
  if (m_Type == "slip") {
    return ("patch");
  }
  if (m_Type == "inlet") {
    return ("patch");
  }
  if (m_Type == "outlet") {
    return ("patch");
  }
  return ("patch");
}
