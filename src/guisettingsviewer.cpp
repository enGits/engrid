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
#include <QtGui>

#include "guisettingsviewer.h"
#include "guisettingstab.h"

#include <QTextStream>
#include <stdio.h>

/* Here is how we we get QTextStreams that look like iostreams */
QTextStream cin(stdin, QIODevice::ReadOnly);
QTextStream cout(stdout, QIODevice::WriteOnly);
QTextStream cerr(stderr, QIODevice::WriteOnly);

GuiSettingsViewer::GuiSettingsViewer(QSettings* Set, QWidget *parent) : QDialog(parent)
{
  m_settings = Set;
  organization = Set->organizationName();
  application = Set->applicationName();
//   CreateViewer();
}

GuiSettingsViewer::GuiSettingsViewer(QString org, QString app, QWidget *parent) : QDialog(parent)
{
  organization = org;
  application = app;
  m_settings = new QSettings(org, app);
//   CreateViewer();
}

void GuiSettingsViewer::CreateViewer()
{
  closeButton = new QPushButton(tr("Close"));
  saveButton = new QPushButton(tr("Save"));

  connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));

  QHBoxLayout *bottomLayout = new QHBoxLayout;
  bottomLayout->addStretch();
  bottomLayout->addWidget(saveButton);
  bottomLayout->addWidget(closeButton);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(&tabWidget);
  mainLayout->addLayout(bottomLayout);
  setLayout(mainLayout);

  setWindowTitle(tr("Settings Viewer"));
  readSettings();
}

void GuiSettingsViewer::save()
{
  for (int i = 0; i < tabWidget.count(); i++) {
    QString group = tabWidget.tabText(i);
    GuiSettingsTab* ST = (GuiSettingsTab*)(tabWidget.widget(i));

    int N;
    QString key;

    if (group != QObject::tr("General")) m_settings->beginGroup(group);

    N = (ST->spinbox_name).size();
    for (int i = 0; i < N; i++) {
      m_settings->beginGroup(QObject::tr("int"));
      key = ST->spinbox_name[i];
      int value = ST->spinbox[i]->value();
      m_settings->setValue(key, value);
      m_settings->endGroup();
    }

    N = (ST->checkbox_name).size();
    for (int i = 0; i < N; i++) {
      m_settings->beginGroup(QObject::tr("bool"));
      key = ST->checkbox_name[i];
      Qt::CheckState value = ST->checkbox[i]->checkState();
      m_settings->setValue(key, value);
      m_settings->endGroup();
    }

    N = (ST->double_lineedit_name).size();
    for (int i = 0; i < N; i++) {
      m_settings->beginGroup(QObject::tr("double"));
      key = ST->double_lineedit_name[i];
      double value = (ST->double_lineedit[i]->text()).toDouble();
      m_settings->setValue(key, value);
      m_settings->endGroup();
    }

    N = (ST->string_lineedit_name).size();
    for (int i = 0; i < N; i++) {
      m_settings->beginGroup(QObject::tr("QString"));
      key = ST->string_lineedit_name[i];
      QString value = ST->string_lineedit[i]->text();
      m_settings->setValue(key, value);
      m_settings->endGroup();
    }

//     N = (ST->filename_dialoglineedit_name).size();
//     for (int i = 0; i < N; i++) {
//       m_settings->beginGroup(QObject::tr("Filename"));
//       key = ST->filename_dialoglineedit_name[i];
//       QString value = ST->filename_dialoglineedit[i]->text();
//       m_settings->setValue(key, value);
//       m_settings->endGroup();
//     }

//     N = (ST->directory_dialoglineedit_name).size();
//     for (int i = 0; i < N; i++) {
//       m_settings->beginGroup(QObject::tr("Directory"));
//       key = ST->directory_dialoglineedit_name[i];
//       QString value = ST->directory_dialoglineedit[i]->text();
//       m_settings->setValue(key, value);
//       m_settings->endGroup();
//     }

    if (group != QObject::tr("General")) m_settings->endGroup();
  }

  close();

}

void GuiSettingsViewer::open()
{
  QDialog dialog(this);

  QLabel *orgLabel = new QLabel(tr("&Organization:"));
  QLineEdit *orgLineEdit = new QLineEdit(organization);
  orgLabel->setBuddy(orgLineEdit);

  QLabel *appLabel = new QLabel(tr("&Application:"));
  QLineEdit *appLineEdit = new QLineEdit(application);
  appLabel->setBuddy(appLineEdit);

  QPushButton *okButton = new QPushButton(tr("OK"));
  okButton->setDefault(true);
  QPushButton *cancelButton = new QPushButton(tr("Cancel"));

  connect(okButton, SIGNAL(clicked()), &dialog, SLOT(accept()));
  connect(cancelButton, SIGNAL(clicked()), &dialog, SLOT(reject()));

  QHBoxLayout *buttonLayout = new QHBoxLayout;
  buttonLayout->addStretch();
  buttonLayout->addWidget(okButton);
  buttonLayout->addWidget(cancelButton);

  QGridLayout *gridLayout = new QGridLayout;
  gridLayout->addWidget(orgLabel, 0, 0);
  gridLayout->addWidget(orgLineEdit, 0, 1);
  gridLayout->addWidget(appLabel, 1, 0);
  gridLayout->addWidget(appLineEdit, 1, 1);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(gridLayout);
  mainLayout->addLayout(buttonLayout);
  dialog.setLayout(mainLayout);

  dialog.setWindowTitle(tr("Choose Settings"));

  if (dialog.exec()) {
    organization = orgLineEdit->text();
    application = appLineEdit->text();
    readSettings();
  }
}

void GuiSettingsViewer::readSettings()
{
  addChildSettings();
  setWindowTitle(tr("Settings Viewer - %1 by %2").arg(application).arg(organization));
}

void GuiSettingsViewer::addChildSettings()
{
  ///\todo Delete for real
  //This only removes the tabs, but does not delete them!!!
  tabWidget.clear();

  tabWidget.addTab(new GuiSettingsTab(organization, application, tr("General")), tr("General"));
  foreach(QString group, m_settings->childGroups()) {
    if ((group != tr("int"))
        && (group != tr("bool"))
        && (group != tr("double"))
        && (group != tr("QString"))
        && (group != tr("Filename"))
        && (group != tr("Directory"))) {
      tabWidget.addTab(new GuiSettingsTab(organization, application, group), group);
    };
  }
}

int GuiSettingsViewer::getSet(QString group, QString key, int value)
{
  QString typed_key = tr("int/") + key;
  if (group != tr("General")) m_settings->beginGroup(group);
  //if key=value pair not found in m_settings file, write it
  if (!m_settings->contains(typed_key)) m_settings->setValue(typed_key, value);
  //read key value from m_settings file and assign it to variable
  int variable = (m_settings->value(typed_key)).toInt();
  if (group != tr("General")) m_settings->endGroup();
  return(variable);
}

double GuiSettingsViewer::getSet(QString group, QString key, double value)
{
  QString typed_key = tr("double/") + key;
  if (group != tr("General")) m_settings->beginGroup(group);
  //if key=value pair not found in m_settings file, write it
  if (!m_settings->contains(typed_key)) m_settings->setValue(typed_key, value);
  //read key value from m_settings file and assign it to variable
  double variable = (m_settings->value(typed_key)).toDouble();
  if (group != tr("General")) m_settings->endGroup();
  return(variable);
}

bool GuiSettingsViewer::getSet(QString group, QString key, bool value)
{
  QString typed_key = tr("bool/") + key;
  if (group != tr("General")) m_settings->beginGroup(group);
  Qt::CheckState state = (Qt::CheckState)(value ? 2 : 0);
  //if key=value pair not found in m_settings file, write it
  if (!m_settings->contains(typed_key)) m_settings->setValue(typed_key, state);
  //read key value from m_settings file and assign it to variable
  bool variable = (m_settings->value(typed_key)).toBool();
  if (group != tr("General")) m_settings->endGroup();
  return(variable);
}

QString GuiSettingsViewer::getSet(QString group, QString key, QString value, int type)
{
  QString typed_key;
  if (type == 0) {
    typed_key = tr("QString/") + key;
  }
  else if (type == 1) {
    typed_key = tr("Filename/") + key;
  }
  else {
    typed_key = tr("Directory/") + key;
  }
  if (group != tr("General")) m_settings->beginGroup(group);
  //if key=value pair not found in m_settings file, write it
  if (!m_settings->contains(typed_key)) m_settings->setValue(typed_key, value);
  //read key value from m_settings file and assign it to variable
  QString variable = (m_settings->value(typed_key)).toString();
  if (group != tr("General")) m_settings->endGroup();
  return(variable);
}
