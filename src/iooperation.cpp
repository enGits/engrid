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
#include "iooperation.h"
#include "guimainwindow.h"

#include <QFileDialog>

IOOperation::IOOperation()
{
  EG_TYPENAME;
  setResetOperationCounter(true);
  setQuickSave(true);
}


void IOOperation::setFormat(QString format)
{
  format_txt = format;  
}

void IOOperation::setExtension(QString extension)
{
  extension_txt = extension;  
}

void IOOperation::readInputFileName()
{
  filename = QFileDialog::getOpenFileName(NULL,"read file",GuiMainWindow::getCwd(),format_txt);
  if (!filename.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(filename).absolutePath());
    valid = true;
  } else {
    valid = false;
  }
}

void IOOperation::readOutputFileName()
{
  filename = QFileDialog::getSaveFileName(NULL,"write file",GuiMainWindow::getCwd(),format_txt);
  if (!filename.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(filename).absolutePath());
    if (filename.right(4) != extension_txt.toLower()) {
      if (filename.right(4) != extension_txt.toUpper()) {
        filename += extension_txt.toLower();
      }
    }
    valid = true;
  } else {
    valid = false;
  }
}

void IOOperation::readOutputDirectory()
{
  filename = QFileDialog::getExistingDirectory(NULL, "write OpenFOAM mesh",GuiMainWindow::getCwd());
  if (!filename.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(filename).absolutePath());
    valid = true;
  } else {
    valid = false;
  }
}

bool IOOperation::isValid()
{
  return valid;
}

char* IOOperation::getCFileName()
{
  return filename.toAscii().data();
}

QString IOOperation::getFileName()
{
  return filename;
}


