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
#include "iooperation.h"
#include "guimainwindow.h"

#include <QFileDialog>

IOOperation::IOOperation()
{
  EG_TYPENAME;
  setResetOperationCounter(true);
  setQuickSave(true);
  m_FileName = "";
}

void IOOperation::setFormat(QString format)
{
  m_FormatTxt = format;
}

void IOOperation::setExtension(QString extension)
{
  m_ExtensionTxt = extension;
}

void IOOperation::readInputFileName(QString default_filename, bool reset)
{
  QApplication::restoreOverrideCursor();
  
  QFileDialog dialog(NULL, "read file", GuiMainWindow::getCwd(), m_FormatTxt);
  dialog.selectFile(default_filename);
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_FileName = selected_files[0];
    if (!m_FileName.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_FileName).absolutePath());
      GuiMainWindow::setUnsaved(true);
      GuiMainWindow::pointer()->setFilename(m_FileName);
      if (reset) {
        GuiMainWindow::pointer()->resetXmlDoc();
      }
      m_Valid = true;
    } else {
      m_Valid = false;
    }
  }
  else m_Valid = false;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void IOOperation::readOutputFileName(QString default_filename)
{
  QApplication::restoreOverrideCursor();
  
  QFileDialog dialog(NULL, "write file", GuiMainWindow::getCwd(), m_FormatTxt);
  dialog.selectFile(default_filename);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);
  if (dialog.exec()) {
    QStringList selected_files = dialog.selectedFiles();
    m_FileName = selected_files[0];
    if (!m_FileName.isNull()) {
      GuiMainWindow::setCwd(QFileInfo(m_FileName).absolutePath());
      if (m_FileName.right(4) != m_ExtensionTxt.toLower()) {
        if (m_FileName.right(4) != m_ExtensionTxt.toUpper()) {
          m_FileName += m_ExtensionTxt.toLower();
        }
      }
      m_Valid = true;
    } else {
      m_Valid = false;
    }
  }
  else m_Valid = false;

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void IOOperation::readOutputDirectory()
{
  QApplication::restoreOverrideCursor();
  m_FileName = QFileDialog::getExistingDirectory(NULL, "write OpenFOAM mesh", GuiMainWindow::getCwd());
  if (!m_FileName.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(m_FileName).absolutePath());
    m_Valid = true;
  } else {
    m_Valid = false;
  }
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

void IOOperation::readInputDirectory(QString title_txt)
{
  QApplication::restoreOverrideCursor();
  m_FileName = QFileDialog::getExistingDirectory(NULL, title_txt, GuiMainWindow::getCwd());
  if (!m_FileName.isNull()) {
    GuiMainWindow::setCwd(QFileInfo(m_FileName).absolutePath());
    m_Valid = true;
  } else {
    m_Valid = false;
  }
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
}

bool IOOperation::isValid()
{
  return m_Valid;
}

const char* IOOperation::getCFileName()
{
  //NOTE: The char pointer will be invalid after the statement in which qPrintable() is used.
  qstrcpy(m_FileName_cc,qPrintable(m_FileName));
  return const_cast<const char *>(m_FileName_cc); //type cast needed for proper output. Technically it will be constant until the next change...
}

QString IOOperation::getFileName()
{
  return m_FileName;
}
