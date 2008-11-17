//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#ifndef iooperation_H
#define iooperation_H

class IOOperation;

#include "operation.h"

#include <QString>
#include <QFile>
#include <QTextStream>

#include <vtkLongArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vector>

class IOOperation : public Operation
{
  
private: // attributes
  
  /** flag to determine if a valid file has been selected */
  bool valid;
  
  /** the file name to read -- normally set by inputReadFileName() */
  QString filename;
  
  /** the file format string (e.g. *.stl, *.vtu, ...) */
  QString format_txt;
  
  /** the file extension for write operations */
  QString extension_txt;

protected: // methods
  
  /**
   * Set the file format string for this Reader.
   * @param format the file name extension (e.g. *.stl, *.vtu, ...)
   */
  void setFormat(QString format);
  
  /**
   * Set the file extension string for write operations.
   * @param extension the file name extension (e.g. *.stl, *.vtu, ...)
   */
  void setExtension(QString format);
  
  /**
   * Get a standard C string representing the file name.
   * @return the file name
   */
  char* getCFileName();
  
  /**
   * Get the file name.
   * @return the file name (as QString)
   */
  QString getFileName();
  
  /**
   * Access to the valid flag (true if a valid file name has been selected)
   * @return the valid flag
   */
  bool isValid();
  
public: // methods
  
  virtual ~IOOperation() {};
  
  /// Open a QFileDialog and make the user input a file name for opening or importing. */
  void readInputFileName();
  
  /// Open a QFileDialog and make the user input a file name for saving or exporting. */
  void readOutputFileName();
  
  /// Open a QFileDialog and make the user input a directory name for exporting. */
  void readOutputDirectory();
  
};

#endif
