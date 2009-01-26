#include "guisettingstab.h"
#include <typeinfo>
#include <QFormLayout>

#include <iostream>
using namespace std;

SettingsTab::SettingsTab(QString org,QString app,QString group,QWidget *parent ): QWidget(parent)
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
