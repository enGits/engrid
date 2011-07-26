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
#include "dialoglineedit.h"

#include <QtGui>
#include <QHBoxLayout>
#include <QtDebug>
#include <QFileInfo>

#include "utilities.h"

DialogLineEdit::DialogLineEdit(QWidget *parent)
    : QWidget(parent)
{
  build();
}

DialogLineEdit::DialogLineEdit(bool openfile, QWidget *parent)
    : QWidget(parent)
{
  build();
  m_openfile = openfile;
}

void DialogLineEdit::build()
{
  setWindowTitle(tr("lineedit with file/directory selection dialog"));
//     this->setHeight(25);
  QHBoxLayout* layout = new QHBoxLayout();
  layout->setSpacing(-1);
  this->setContentsMargins(0, 0, 0, 0);
  layout->setContentsMargins(0, 0, 0, 0);
  this->setLayout(layout);
  m_button = new QPushButton(tr("..."));
  m_lineedit = new QLineEdit();
  m_button->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
  m_lineedit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
  layout->addWidget(m_lineedit);
  layout->addWidget(m_button);

  connect(m_button, SIGNAL(clicked()), this, SLOT(openDialog()));

  m_button->setMaximumWidth(16);
  this->resize(180, 25);
//    this->resize(0, 0);
// setGeometry(10, 40, 180, 40);

  // by default, open directory and use file dialog
  m_openfile = false;
  m_UseFileDialog = true;
}

void DialogLineEdit::openDialog()
{
  emit clicked();

  if (m_UseFileDialog) {
    QString name;
//     qDebug() << "m_dir=" << m_dir;
    if (m_openfile) name = QFileDialog::getOpenFileName(NULL, m_caption, m_dir, m_filter);
    else name = getDirectory(NULL, m_caption, m_dir);

    if (!name.isNull()) {
      m_lineedit->setText(name);
      m_dir = name;//getDir(name);
    }
  }
}

QString DialogLineEdit::text()
{
  return(m_lineedit->text());
}

void DialogLineEdit::setText(QString str)
{
  m_lineedit->setText(str);
  m_dir = str;//getDir(str);
}

bool DialogLineEdit::getMode()
{
  return m_openfile;
}

void DialogLineEdit::setMode(bool mode)
{
  m_openfile = mode;
}

void DialogLineEdit::setCaption(QString str)
{
  m_caption = str;
}

void DialogLineEdit::setDir(QString str)
{
  m_dir = str;
}

void DialogLineEdit::setFilter(QString str)
{
  m_filter = str;
}

QString DialogLineEdit::getCaption()
{
  return m_caption;
}

QString DialogLineEdit::getDir()
{
  return m_dir;
}

QString DialogLineEdit::getFilter()
{
  return m_filter;
}

// QString DialogLineEdit::getDir(QString str)
// {
//   if (m_openfile) {
//     return str;
//   }
//   else {
//     QFileInfo fileinfo(str);
//     return fileinfo.absolutePath();
//   }
// }

void DialogLineEdit::useFileDialog(bool arg)
{
  m_UseFileDialog = arg;
}
