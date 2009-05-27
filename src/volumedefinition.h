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
#ifndef VOLUMEDEFINITION_H
#define VOLUMEDEFINITION_H

#include <QString>
#include <QMap>

class VolumeDefinition
{

  QString name;
  int vc;
  QMap<int,int> bcs;

public: // methods

  VolumeDefinition();
  VolumeDefinition(QString name, int vc);

  QString getName() { return name; }
  int getVC() { return vc; }
  void addBC(int bc, int sign);
  int getSign(int bc);

};

#endif // VOLUMEDEFINITION_H
