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
#include "guisettingstab.h"
#include <typeinfo>
#include <QFormLayout>
#include <QtDebug>

#include <iostream>
using namespace std;

GuiSettingsTab::GuiSettingsTab(QString org, QString app, QString group, QWidget *parent): QWidget(parent)
{

//   qDebug() << "group=" << group;

  QFormLayout *permissionsLayout = new QFormLayout;

  QSettings settings(org, app);

  if (group != QObject::tr("General")) settings.beginGroup(group);

  settings.beginGroup(QObject::tr("int"));
  foreach(QString key, settings.childKeys()) {
    int I = settings.value(key).toInt();
    spinbox.append(new QSpinBox);
    spinbox_name.append(key);
    spinbox.back()->setValue(I);
    permissionsLayout->addRow(key, spinbox.back());
  }
  settings.endGroup();

  settings.beginGroup(QObject::tr("bool"));
  foreach(QString key, settings.childKeys()) {
    checkbox.append(new QCheckBox);
    checkbox_name.append(key);
    checkbox.back()->setCheckState((Qt::CheckState)settings.value(key).toInt());
    permissionsLayout->addRow(key, checkbox.back());
  }
  settings.endGroup();

  settings.beginGroup(QObject::tr("double"));
  foreach(QString key, settings.childKeys()) {
    double D = settings.value(key).toDouble();
    QString s;
    s.setNum(D);
    //s.sprintf("%.2f", D);
    double_lineedit.append(new QLineEdit);
    double_lineedit_name.append(key);
    double_lineedit.back()->setText(s);
    permissionsLayout->addRow(key, double_lineedit.back());
  }
  settings.endGroup();

  settings.beginGroup(QObject::tr("QString"));
  foreach(QString key, settings.childKeys()) {
    QString S = settings.value(key).toString();
    string_lineedit.append(new QLineEdit);
    string_lineedit_name.append(key);
    string_lineedit.back()->setText(S);
    permissionsLayout->addRow(key, string_lineedit.back());
  }
  settings.endGroup();

  settings.beginGroup(QObject::tr("Filename"));
  foreach(QString key, settings.childKeys()) {
    QString S = settings.value(key).toString();
    filename_dialoglineedit.append(new DialogLineEdit(true));
    filename_dialoglineedit_name.append(key);
    filename_dialoglineedit.back()->setText(S);
    permissionsLayout->addRow(key, filename_dialoglineedit.back());
  }
  settings.endGroup();

  settings.beginGroup(QObject::tr("Directory"));
  foreach(QString key, settings.childKeys()) {
    QString S = settings.value(key).toString();
    directory_dialoglineedit.append(new DialogLineEdit(false));
    directory_dialoglineedit_name.append(key);
    directory_dialoglineedit.back()->setText(S);
    permissionsLayout->addRow(key, directory_dialoglineedit.back());
  }
  settings.endGroup();

  if (group != QObject::tr("General")) settings.endGroup();

  setLayout(permissionsLayout);

}
