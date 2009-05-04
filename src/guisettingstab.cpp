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
#include "guisettingstab.h"
#include <typeinfo>
#include <QFormLayout>

#include <iostream>
using namespace std;

GuiSettingsTab::GuiSettingsTab(QString org, QString app, QString group, QWidget *parent): QWidget(parent)
{
  
  QFormLayout *permissionsLayout = new QFormLayout;

  QSettings settings(org, app);
  
  if(group!="General") settings.beginGroup(group);
  
  settings.beginGroup("int");
  foreach (QString key, settings.childKeys()) {
    int I=settings.value(key).toInt();
    spinbox.append(new QSpinBox);
    spinbox_name.append(key);
    spinbox.back()->setValue(I);
    permissionsLayout->addRow(key, spinbox.back());
  }
  settings.endGroup();
  
  settings.beginGroup("bool");
  foreach (QString key, settings.childKeys()) {
    checkbox.append(new QCheckBox);
    checkbox_name.append(key);
    checkbox.back()->setCheckState((Qt::CheckState)settings.value(key).toInt());
    permissionsLayout->addRow(key, checkbox.back());
  }
  settings.endGroup();
  
  settings.beginGroup("double");
  foreach (QString key, settings.childKeys()) {
    double D=settings.value(key).toDouble();
    QString s;
    s.sprintf("%.2f", D);
    lineedit.append(new QLineEdit);
    lineedit_name.append(key);
    lineedit.back()->setText(s);
    permissionsLayout->addRow(key, lineedit.back());
  }
  settings.endGroup();

  if(group!="General") settings.endGroup();
  
  setLayout(permissionsLayout);
  
}
