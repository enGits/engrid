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
#include "physicalboundaryconditions.h"

#include <QStringList>
#include <QtDebug>

PhysicalBoundaryConditions::PhysicalBoundaryConditions()
{
  this->m_Name = "unknown";
  this->m_Index = -1;
  setDefaults();
}

PhysicalBoundaryConditions::PhysicalBoundaryConditions(QString name, int index)
{
  this->m_Name = name;
  this->m_Index = index;
  setDefaults();
}

PhysicalBoundaryConditions::PhysicalBoundaryConditions(QString name, int index, QString values)
{
  this->m_Name = name;
  this->m_Index = index;
  setDefaults();
  
  qWarning()<<"values="<<values;
  
  QStringList L = values.split(";");
  qWarning()<<"L="<<L;
  for(int i=0;i<L.size();i++) {
    QStringList L_pair = L[i].split("=");
    if(L_pair[0].trimmed()=="Pressure") m_Pressure = L_pair[1].toDouble();
    if(L_pair[0].trimmed()=="Temperature") m_Temperature = L_pair[1].toDouble();
    if(L_pair[0].trimmed()=="Velocity") m_Velocity = L_pair[1].toDouble();
  }
}

void PhysicalBoundaryConditions::setDefaults()
{
  m_Pressure = 1;
  m_Temperature = 2;
  m_Velocity = 3;
}

QString PhysicalBoundaryConditions::getIndex()
{
  QString ret;
  ret.setNum(m_Index);
  return(ret);
}

QString PhysicalBoundaryConditions::getName()
{
  return(m_Name);
}

QString PhysicalBoundaryConditions::getValues()
{
  QString ret("");
  QString str_Pressure; str_Pressure.setNum(m_Pressure); ret += "Pressure=" + str_Pressure + ";";
  QString str_Temperature; str_Temperature.setNum(m_Temperature); ret += "Temperature=" + str_Temperature + ";";
  QString str_Velocity; str_Velocity.setNum(m_Velocity); ret += "Velocity=" + str_Velocity + ";";
  return(ret);
}
