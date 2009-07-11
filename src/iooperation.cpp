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
  m_filename = "";
}

void IOOperation::setFormat(QString format)
{
  m_format_txt = format;
}

void IOOperation::setExtension(QString extension)
{
  m_extension_txt = extension;
}

///@@@ TODO: Low priority: add default name argument so that it doesn't offer a file/directory of the wrong type by default. Must be done in GuiMainWindow too.
void IOOperation::readInputFileName()
{
  QFileDialog dialog(NULL,"read file",GuiMainWindow::getCwd(),m_format_txt);
  dialog.selectFile(GuiMainWindow::pointer()->getFilename());
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_filename = selected_files[0];
    if (!m_filename.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_filename).absolutePath());
      GuiMainWindow::setUnsaved(true);
      GuiMainWindow::pointer()->setFilename(m_filename);
      m_valid = true;
    } else {
      m_valid = false;
    }
  }
  else m_valid = false;
}

void IOOperation::readOutputFileName()
{
  QFileDialog dialog(NULL,"write file",GuiMainWindow::getCwd(),m_format_txt);
  dialog.selectFile(GuiMainWindow::pointer()->getFilename());
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_filename = selected_files[0];
    if (!m_filename.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_filename).absolutePath());
      if (m_filename.right(4) != m_extension_txt.toLower()) {
        if (m_filename.right(4) != m_extension_txt.toUpper()) {
          m_filename += m_extension_txt.toLower();
        }
      }
      m_valid = true;
    } else {
      m_valid = false;
    }
  }
  else m_valid = false;
}

void IOOperation::readOutputDirectory()
{
  QFileDialog dialog(NULL, "write OpenFOAM mesh",GuiMainWindow::getCwd());
  dialog.selectFile(GuiMainWindow::pointer()->getFilename());
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_filename = selected_files[0];
    if (!m_filename.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_filename).absolutePath());
      m_valid = true;
    } else {
      m_valid = false;
    }
  }
  else m_valid = false;
}

void IOOperation::readInputDirectory(QString title_txt)
{
  QFileDialog dialog(NULL, title_txt, GuiMainWindow::getCwd());
  dialog.selectFile(GuiMainWindow::pointer()->getFilename());
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_filename = selected_files[0];
    if (!m_filename.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_filename).absolutePath());
      GuiMainWindow::setUnsaved(true);
      GuiMainWindow::pointer()->setFilename(m_filename);
      m_valid = true;
    } else {
      m_valid = false;
    }
  }
  else m_valid = false;
}

bool IOOperation::isValid()
{
  return m_valid;
}

const char* IOOperation::getCFileName()
{
  return qPrintable(m_filename);
}

QString IOOperation::getFileName()
{
  return m_filename;
}
