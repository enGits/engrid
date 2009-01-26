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
