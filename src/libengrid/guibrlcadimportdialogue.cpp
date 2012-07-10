//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "guibrlcadimportdialogue.h"
#include "ui_guibrlcadimportdialogue.h"

#include <QProcess>

GuiBrlCadImportDialogue::GuiBrlCadImportDialogue(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::GuiBrlCadImportDialogue)
{
  ui->setupUi(this);
}

GuiBrlCadImportDialogue::~GuiBrlCadImportDialogue()
{
  delete ui;
}

void GuiBrlCadImportDialogue::prepare(QString file_name)
{
  QString program = "/usr/brlcad/bin/mged";
  QStringList arguments;
  QProcess proc(this);
  arguments << "-c" << file_name<< "ls";
  proc.start(program, arguments);
  proc.waitForFinished();
  QString output = proc.readAllStandardOutput() + proc.readAllStandardError();
  QStringList objects = output.split(QRegExp("\\s+"));
  ui->listWidget->clear();
  foreach (QString obj, objects) {
    ui->listWidget->addItem(obj);
  }
}

bool GuiBrlCadImportDialogue::hasSelectedObject()
{
  if (ui->listWidget->selectedItems().size() == 1) {
    return true;
  }
  return false;
}

QString GuiBrlCadImportDialogue::selectedObject()
{
  if (hasSelectedObject()) {
    QString object_txt = ui->listWidget->selectedItems().first()->text();
    QString object = object_txt.split("/").first();
    return object;
  } else {
    EG_BUG;
  }
}

double GuiBrlCadImportDialogue::scanMemory()
{
  double mem = 1024.0;
  mem *= mem*mem;
  mem *= ui->doubleSpinBoxMemory->value();
  return mem;
}



