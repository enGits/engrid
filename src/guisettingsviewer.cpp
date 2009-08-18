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
  settings = Set;
  organization = Set->organizationName();
  application = Set->applicationName();
//   CreateViewer();
}

GuiSettingsViewer::GuiSettingsViewer(QString org, QString app, QWidget *parent) : QDialog(parent)
{
  organization = org;
  application = app;
  settings = new QSettings(org, app);
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

    if (group != "General") settings->beginGroup(group);

    N = (ST->spinbox_name).size();
    for (int i = 0; i < N; i++) {
      settings->beginGroup("int");
      key = ST->spinbox_name[i];
      int value = ST->spinbox[i]->value();
      settings->setValue(key, value);
      settings->endGroup();
    }

    N = (ST->checkbox_name).size();
    for (int i = 0; i < N; i++) {
      settings->beginGroup("bool");
      key = ST->checkbox_name[i];
      Qt::CheckState value = ST->checkbox[i]->checkState();
      settings->setValue(key, value);
      settings->endGroup();
    }

    N = (ST->lineedit_name).size();
    for (int i = 0; i < N; i++) {
      settings->beginGroup("double");
      key = ST->lineedit_name[i];
      double value = (ST->lineedit[i]->text()).toDouble();
      settings->setValue(key, value);
      settings->endGroup();
    }

    if (group != "General") settings->endGroup();
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
  ///@@@  TODO: Delete for real
  //This only removes the tabs, but does not delete them!!!
  tabWidget.clear();

  tabWidget.addTab(new GuiSettingsTab(organization, application, "General"), "General");
  foreach(QString group, settings->childGroups()) {
    if ((group != "int") && (group != "bool") && (group != "double")) {
      tabWidget.addTab(new GuiSettingsTab(organization, application, group), group);
    };
  }
}
