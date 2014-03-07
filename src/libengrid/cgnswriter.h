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
#ifndef CGNSWRITER_H
#define CGNSWRITER_H

class CgnsWriter;

#include "iooperation.h"

#ifdef CGNS_SUPPORT
#include "cgnslib.h"
#endif

/**
 * Writer for CGNS files.
 */
class CgnsWriter : public IOOperation
{

protected: // attributes

  int fn;
  int B;
  int Z;
  QVector<int> eg2cgns;

protected: // methods

  void writeGrid();
  void writeBcs();
  virtual void operate();

public: // methods

  CgnsWriter();

};

#endif // CGNSWRITER_H
