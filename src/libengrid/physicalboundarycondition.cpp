// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "physicalboundarycondition.h"

#include <QStringList>
#include <QtDebug>

PhysicalBoundaryCondition::PhysicalBoundaryCondition()
{
  m_Name = "unknown";
  m_Index = -1;
  m_Type = "unknown";
}

void PhysicalBoundaryCondition::setType(QString type)
{
  m_Type = type;
  m_VarNames.clear();
  m_VarValues.clear();
  if (m_Type == "symmetry") {
  }
  if (m_Type == "wall") {
  }
  if (m_Type == "slip") {
  }
  if (m_Type == "inlet") {
    m_VarNames.push_back("velocity");
    m_VarValues.push_back(0);
    m_VarNames.push_back("turbulent-intensity");
    m_VarValues.push_back(0.04);
    m_VarNames.push_back("turbulent-length-scale");
    m_VarValues.push_back(1);
    m_VarNames.push_back("temperature");
    m_VarValues.push_back(300);
  }
  if (m_Type == "outlet") {
    m_VarNames.push_back("pressure");
    m_VarValues.push_back(0);
  }
}

QString PhysicalBoundaryCondition::getFoamEpsilon()
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
    double k       = 1.5*sqr(getVarValue(0)*getVarValue(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValue(2);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << epsilon << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamK()
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
    double k       = 1.5*sqr(getVarValue(0)*getVarValue(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValue(2);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << k << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamOmega()
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
    double k       = 1.5*sqr(getVarValue(0)*getVarValue(1));
    double epsilon = (pow(0.09, 0.75)*pow(k, 1.5))/getVarValue(2);
    double omega   = epsilon/(0.09*k);
    s << "        type  fixedValue;\n";
    s << "        value uniform " << omega << ";\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamP()
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
    s << "        value uniform " << getVarValue(0) << ";\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamU(vec3_t n)
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
    s << "        value uniform (" << getVarValue(0)*n[0] << " " << getVarValue(0)*n[1] << " " << getVarValue(0)*n[2] << ");\n";
  }
  if (m_Type == "outlet") {
    s << "        type zeroGradient;\n";
  }
  return str;
}

QString PhysicalBoundaryCondition::getFoamT()
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
    s << "        value uniform " << getVarValue(3) << ";\n";
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
