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
#ifndef SETTINGSTAB_H
#define SETTINGSTAB_H

#include <QWidget>
#include <QString>
#include <QtGui>
#include <QVector>

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class SettingsTab : public QWidget
{
Q_OBJECT
public:
     QVector<QSpinBox*> spinbox;

public:
	//constructors
    SettingsTab(QWidget *parent = 0);
    SettingsTab(QString org,QString app,QString group,QWidget *parent = 0);
};

#endif
