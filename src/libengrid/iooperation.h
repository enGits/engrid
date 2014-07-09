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
#ifndef iooperation_H
#define iooperation_H

class IOOperation;

#include "operation.h"
// #include "guimainwindow.h"

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QFileInfo>

#include <vtkLongArray.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vector>

class IOOperation : public Operation
{
  
private: // attributes
  
  bool    m_Valid;          ///< flag to determine if a valid file has been selected
  QString m_FileName;       ///< file name to read -- normally set by inputReadFileName()
  char m_FileName_cc[1024]; ///< "const char *" copy version of file name to read -- only used by getCFileName()
  QString m_FormatTxt;      ///< file format string (e.g. *.stl, *.vtu, ...)
  QString m_ExtensionTxt;   ///< file extension for write operations
  bool    m_FileNameSet;    ///< has the file name been set already (batch operation)?

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
  const char* getCFileName();
  
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
  
  IOOperation();
  virtual ~IOOperation() {}

  void readInputFileName(QString default_filename, bool reset = true);   ///< Open a QFileDialog and make the user input a file name for opening or importing.
  void readOutputFileName(QString default_filename);                     ///< Open a QFileDialog and make the user input a file name for saving or exporting.
  void readOutputDirectory();                                            ///< Open a QFileDialog and make the user input a directory name for exporting.
  void readInputDirectory(QString title_txt = "select input directory"); ///< Open a QFileDialog and make the user input a directory name for importing.
  void setFileName(QString file_name);                                   ///< set the file name for batch operation
  
};

#endif
