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
 * Creates a QWidget listing all key/value pairs contained in the group "group" 
 * of the QSettings file corresponding to the (org,app) pair.
 * integers appear in spinboxes
 * doubles appear in line edit boxes
 * booleans appear in checkboxes
 */
class GuiSettingsTab : public QWidget
{
  
  Q_OBJECT;
    
public:
  
  QVector<QString> spinbox_name;
  QVector<QSpinBox*> spinbox;
  
  QVector<QString> checkbox_name;
  QVector<QCheckBox*> checkbox;
  
  QVector<QString> lineedit_name;
  QVector<QLineEdit*> lineedit;
  
public:
	//constructors
  /**
   * Constructor using the (org,app) pair to determine QSettings
   * @param org organization
   * @param app application
   * @param group group
   * @param parent Parent QWidget
   */
  GuiSettingsTab(QString org,QString app,QString group,QWidget *parent = 0);
  
};

#endif
