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
#ifndef PHYSICALBOUNDARYCONDITIONS_H
#define PHYSICALBOUNDARYCONDITIONS_H

#include <QString>

class PhysicalBoundaryConditions{
public:
  QString m_Name;
  int m_Index;
public:
  double m_Pressure;
  double m_Temperature;
  double m_Velocity;
public:
  PhysicalBoundaryConditions();
  PhysicalBoundaryConditions(QString name, int index);
  PhysicalBoundaryConditions(QString name, int index, QString values);
public:
  QString getName();
  QString getIndex();
  QString getValues();
  void setDefaults();
};

#endif
