#include "guisettingstab.h"
#include <typeinfo>
#include <QFormLayout>

#include <iostream>
using namespace std;

SettingsTab::SettingsTab(QWidget *parent)
 : QWidget(parent)
{
	cout<<"[SettingsTab::SettingsTab(QWidget *parent)]"<<endl;
}

SettingsTab::SettingsTab(QString org,QString app,QString group,QWidget *parent ): QWidget(parent)
{
  cout<<"[SettingsTab::SettingsTab(QString str, QWidget *parent) : QWidget(parent)]"<<endl;
  
  QFormLayout *permissionsLayout = new QFormLayout;

  QSettings settings(org, app);
  
  if(group!="General") settings.beginGroup(group);
  
  cout<<"++++++++++++++++++++++++++"<<endl;
  cout<<"GROUP="<<group.toLatin1().data()<<endl;
  cout<<"++++++++++++++++++++++++++"<<endl;
  settings.beginGroup("int");
  foreach (QString key, settings.childKeys()) {
    int I=settings.value(key).toInt();
    cout << "int: " << key.toLatin1().data() << "="<< I << endl;
    spinbox.append(new QSpinBox);
    spinbox_name.append(key);
    spinbox.back()->setValue(I);
    permissionsLayout->addRow(key, spinbox.back());
  }
  settings.endGroup();
  
  settings.beginGroup("bool");
  foreach (QString key, settings.childKeys()) {
    cout << "bool: "<< key.toLatin1().data() << "="<< settings.value(key).toBool() << endl;
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
    cout << "double: "<< key.toLatin1().data() << "="<< D << endl;
    lineedit.append(new QLineEdit);
    lineedit_name.append(key);
    lineedit.back()->setText(s);
    permissionsLayout->addRow(key, lineedit.back());
  }
  settings.endGroup();
  cout<<"++++++++++++++++++++++++++"<<endl;
  
//   QStringList keys = settings.childKeys();
  
//   int N=keys.size();
//   cout << "number of variables = " << N <<endl;

  cout<<"=============="<<endl;
  foreach (QString key, settings.childKeys()) {
    cout << key.toLatin1().data() << "="<< settings.value(key).toInt() << endl;
  }
  cout<<"=============="<<endl;
  if(group!="General") settings.endGroup();
  
  setLayout(permissionsLayout);
  
}
