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
#ifndef stlreader_H
#define stlreader_H

class StlReader;

#include "iooperation.h"

/**
 * Reader for ASCII and binary STL files
 */
class StlReader : public IOOperation
{

private: // attributes

  double  m_Tolerance;
  bool    m_FileNameSet;
  QString m_FileName;
  int     m_MaxNumCleanIter;
  
protected: // methods
  
  virtual void operate();
  
public: // methods
  
  /** The constructor sets the file format string. */
  StlReader();

  void setTolerance(double tol) { m_Tolerance = tol; }
  void setFileName(QString file_name);
  void setMaximalCleaningIterations(int N) { m_MaxNumCleanIter = N; }
    
};

#endif
