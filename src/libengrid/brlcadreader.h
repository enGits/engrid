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
#ifndef BRLCADREADER_H
#define BRLCADREADER_H

class BrlcadReader;

#include "iooperation.h"

#include <QMap>

class BrlcadReader : public IOOperation
{

private: // attributes

  QList<vtkUnstructuredGrid*> m_Grids;
  QMap<vtkUnstructuredGrid*, QString> m_BCNames;


protected: // methods

  void processStlFile(QString file_name, bool append_to_list = true);
  void findBoundaryCodes();

  virtual void operate();


public: // methods

  BrlcadReader();

};

#endif // BRLCADREADER_H
