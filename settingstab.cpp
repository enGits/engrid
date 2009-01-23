#include "settingstab.h"
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
  
  QStringList keys = settings.childKeys();
  
  int N=keys.size();
  cout << "number of variables = " << N <<endl;

  cout<<"=============="<<endl;
  foreach (QString key, settings.childKeys()) {
    cout << key.toLatin1().data() << "="<< settings.value(key).toInt() << endl;
    spinbox.append(new QSpinBox);
    spinbox.back()->setValue(settings.value(key).toInt());
    permissionsLayout->addRow(key, spinbox.back());
  }
  cout<<"=============="<<endl;
  if(group!="General") settings.endGroup();
  
  setLayout(permissionsLayout);
  
}
