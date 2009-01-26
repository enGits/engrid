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
#ifndef SETTINGSVIEWER_H
#define SETTINGSVIEWER_H

#include <QDialog>
#include <QtGui>

#include "settingstab.h"

// class QPushButton;
// class QSettings;
// class QTreeWidget;
// class QTreeWidgetItem;

class SettingsViewer : public QDialog
{
    Q_OBJECT

public:
	//constructors
  SettingsViewer(QString org, QString app,QWidget *parent = 0);
  SettingsViewer(QSettings* Set,QWidget *parent = 0);
  
private slots:
    void open();
    void save();
    void CreateViewer();
    void readSettings();
    void addChildSettings();

private:
//     QTreeWidget *treeWidget;
    QPushButton *openButton;
    QPushButton *saveButton;
    QPushButton *closeButton;

    QTabWidget tabWidget;
    QVector<SettingsTab> tabs;  
  
    QString organization;
    QString application;
    QSettings* settings;
};

#endif
