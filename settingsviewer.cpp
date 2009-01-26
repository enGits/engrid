#include <QtGui>

#include "settingsviewer.h"
#include "settingstab.h"

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
  cout<<"[SettingsViewer::SettingsViewer(QString org, QString app,QWidget *parent ) : QDialog(parent)]"<<endl;
  
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
  cout<<"=============="<<endl;
  for (int i=0;i<tabWidget.count();i++)
  {
    QString group=tabWidget.tabText(i);
    cout<<"group="<<group.toLatin1().data()<<endl;
    SettingsTab* ST=(SettingsTab*)(tabWidget.widget(i));
    ST->spinbox;
    ST->spinbox_name;
    int N;
    
    if(group!="General") settings->beginGroup(group);
    
    N=(ST->spinbox_name).size();
    cout<<"N1="<<N<<endl;
    QString key;
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("int");
      cout<<"name="<<ST->spinbox_name[i]<<endl;
      cout<<"value="<<ST->spinbox[i]->value()<<endl;
      key=ST->spinbox_name[i];
      int value=ST->spinbox[i]->value();
      settings->setValue(key,value);
      settings->endGroup();
    }
    N=(ST->checkbox_name).size();
    cout<<"N2="<<N<<endl;
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("bool");
      cout<<"name="<<ST->checkbox_name[i]<<endl;
      cout<<"value="<<ST->checkbox[i]->checkState()<<endl;
      cout<<"ST->spinbox_name[i]="<<ST->spinbox_name[i]<<endl;
      key=ST->checkbox_name[i];
      Qt::CheckState value=ST->checkbox[i]->checkState();
      settings->setValue(key,value);
      settings->endGroup();
    }
    N=(ST->lineedit_name).size();
    cout<<"N3="<<N<<endl;
    for(int i=0;i<N;i++)
    {
      settings->beginGroup("double");
      cout<<"name="<<ST->lineedit_name[i]<<endl;
      cout<<"value="<<ST->lineedit[i]->text()<<endl;
      key=ST->lineedit_name[i];
      double value=(ST->lineedit[i]->text()).toDouble();
      settings->setValue(key,value);
      settings->endGroup();
    }
    
    if(group!="General") settings->endGroup();
  }
  cout<<"=============="<<endl;

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
//     treeWidget->clear();
    addChildSettings();

/*    treeWidget->sortByColumn(0);
    treeWidget->setFocus();*/
    setWindowTitle(tr("Settings Viewer - %1 by %2").arg(application).arg(organization));
}

void SettingsViewer::addChildSettings()
{
  tabWidget.clear(); //This only removes the tabs, but does not delete them!!! TODO: delete for real

  cout<<"=============="<<endl;
  tabWidget.addTab(new SettingsTab(organization, application, "General"), "General");
  foreach (QString group, settings->childGroups()) {
    tabWidget.addTab(new SettingsTab(organization, application, group), group);
  }
  cout<<"=============="<<endl;
}
