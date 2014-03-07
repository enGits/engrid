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
#ifndef gmshwriter_H
#define gmshwriter_H

class GmshWriter;

#include "gmshiooperation.h"

/**
 * Writer for Gmsh files; this Writer currently only supports
 * version 1.0 of the the Gmsh file format.
 */
class GmshWriter : public GmshIOOperation
{
  
  void writeAscii1(vtkUnstructuredGrid *grid);
  void writeAscii2(vtkUnstructuredGrid *grid);
  
protected: // methods
  
  virtual void operate();
  
public: // methods
  
  /** The constructor sets the file format string. */
  GmshWriter();
  
};

#endif
