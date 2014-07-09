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
#ifndef VOLUMEDEFINITION_H
#define VOLUMEDEFINITION_H

#include <QString>
#include <QMap>

class VolumeDefinition
{
private:
  QString name;      ///< name of the volume (unique identifier)
  int vc;            ///< volume code
  QMap<int,int> bcs; ///< maps each BC of the volume to a sign (+1 if the outside is green /-1 if the outside is yellow)

public: // methods

  VolumeDefinition();
  VolumeDefinition(QString name, int vc);

  void setVC(int v) { vc = v; }

  QString getName() { return name; } ///< get name
  int getVC() { return vc; }         ///< get volume code
  void addBC(int bc, int sign);      ///< bcs[bc] = sign;
  int getSign(int bc);               ///< returns bcs[bc] if it exists, otherwise returns 0
};

#endif // VOLUMEDEFINITION_H
