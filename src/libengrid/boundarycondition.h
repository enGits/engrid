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
#ifndef boundarycondition_H
#define boundarycondition_H

class BoundaryCondition;

#include <QString>

class BoundaryCondition
{
  
private: // attributes
  
  int     m_Code;
  QString m_Name;
  QString m_Type;
  

public: // methods
  
  BoundaryCondition();
  BoundaryCondition(QString name, QString type, int code = -1);
  
  int     getCode() { return m_Code; }
  QString getName() { return m_Name; }
  QString getType() { return m_Type; }

  void setCode(int code) { m_Code = code; }
  void setName(QString name) { m_Name = name; }
  void setType(QString type) { m_Type = type; }

  bool operator==(const BoundaryCondition& bc) const { return m_Name == bc.m_Name && m_Type == bc.m_Type; }
  bool operator<(const BoundaryCondition& bc) const;
  
};

#endif
