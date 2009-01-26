#include <QtGui>

#include "guisettingsviewer.h"
#include "guisettingstab.h"

// #include <iostream>
// using namespace std;

#include <QTextStream>
#include <stdio.h>

/* Here is how we we get QTextStreams that look like iostreams */
QTextStream cin(stdin, QIODevice::ReadOnly);
QTextStream cout(stdout, QIODevice::WriteOnly);
QTextStream cerr(stderr, QIODevice::WriteOnly);

SettingsViewer::SettingsViewer(QSettings* Set,QWidget *parent) : QDialog(parent)
{
  settings = Set;
  organization =Set->organizationName();
  application = Set->applicationName();
  CreateViewer();
}

SettingsViewer::SettingsViewer(QString org, QString app,QWidget *parent ) : QDialog(parent)
{
    organization = org;
    application = app;
    settings=new QSettings(org,app);
    CreateViewer();
}

void SettingsViewer::CreateViewer()
{
//     openButton = new QPushButton(tr("&Open..."));
//     openButton->setDefault(true);
  
  closeButton = new QPushButton(tr("Close"));
  saveButton = new QPushButton(tr("Save"));
  
//     connect(openButton, SIGNAL(clicked()), this, SLOT(open()));
  connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));
  connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
  
  QHBoxLayout *bottomLayout = new QHBoxLayout;
  bottomLayout->addStretch();
//   bottomLayout->addWidget(openButton);
  bottomLayout->addWidget(saveButton);
  bottomLayout->addWidget(closeButton);
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(&tabWidget);
  mainLayout->addLayout(bottomLayout);
  setLayout(mainLayout);
  
  setWindowTitle(tr("Settings Viewer"));
  readSettings();
}

void SettingsViewer::save()
{
  for (int i=0;i<tabWidget.count();i++)
  {
    QString group=tabWidget.tabText(i);
    SettingsTab* ST=(SettingsTab*)(tabWidget.widget(i));
    
    int N;
    QString key;
    
    if(group!="General") settings->beginGroup(group);
    
    N=(ST->spinbox_name).size();
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("int");
      key=ST->spinbox_name[i];
      int value=ST->spinbox[i]->value();
      settings->setValue(key,value);
      settings->endGroup();
    }
    
    N=(ST->checkbox_name).size();
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("bool");
      key=ST->checkbox_name[i];
      Qt::CheckState value=ST->checkbox[i]->checkState();
      settings->setValue(key,value);
      settings->endGroup();
    }
    
    N=(ST->lineedit_name).size();
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("double");
      key=ST->lineedit_name[i];
      double value=(ST->lineedit[i]->text()).toDouble();
      settings->setValue(key,value);
      settings->endGroup();
    }
    
    if(group!="General") settings->endGroup();
  }

}

void SettingsViewer::open()
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

void SettingsViewer::readSettings()
{
    addChildSettings();
    setWindowTitle(tr("Settings Viewer - %1 by %2").arg(application).arg(organization));
}

void SettingsViewer::addChildSettings()
{
  tabWidget.clear(); //This only removes the tabs, but does not delete them!!! TODO: delete for real

  tabWidget.addTab(new SettingsTab(organization, application, "General"), "General");
  foreach (QString group, settings->childGroups()) {
    tabWidget.addTab(new SettingsTab(organization, application, group), group);
  }
}
