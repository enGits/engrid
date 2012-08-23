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
#include <QFileInfo>

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
  m_StlFileName = file_name + ".stl";
  if (QFileInfo(m_StlFileName).exists()) {
    ui->checkBoxSTL->setEnabled(true);
    ui->checkBoxSTL->setChecked(true);
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

int GuiBrlCadImportDialogue::smoothingIterations()
{
  return ui->spinBoxSmoothing->value();
}

int GuiBrlCadImportDialogue::preservationType()
{
  if (ui->radioButtonNoPreservation->isChecked()) return 0;
  if (ui->radioButtonSolidPreservation->isChecked()) return 1;
  return 2;
}

double GuiBrlCadImportDialogue::smallestFeatureSize()
{
  return ui->lineEditSmallestFeature->text().toDouble();
}

double GuiBrlCadImportDialogue::smallestResolution()
{
  return ui->lineEditMinResolution->text().toDouble();
}

bool GuiBrlCadImportDialogue::useStlFile()
{
  return ui->checkBoxSTL->isChecked();
}

QString GuiBrlCadImportDialogue::stlFileName()
{
  return m_StlFileName;
}

double GuiBrlCadImportDialogue::reduction()
{
  return 0.01*double(ui->horizontalSliderReduction->value());
}
