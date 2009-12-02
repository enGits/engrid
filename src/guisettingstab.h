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

#include "dialoglineedit/dialoglineedit.h"
///\todo fix problems with more than one boolean + general tab stuff

/**
 * Creates a QWidget listing all key/value pairs contained in the group "group"
 * of the QSettings file corresponding to the (org,app) pair.
 * integers appear in spinboxes
 * doubles appear in line edit boxes
 * booleans appear in checkboxes
 * \todo fix problems with more than one boolean + general tab stuff
 */
class GuiSettingsTab : public QWidget
{

    Q_OBJECT;

  public:

    QVector<QString> spinbox_name;///< vector of the label texts of the QSpinBox widgets
    QVector<QSpinBox*> spinbox;///< vector of pointers to the QSpinBox widgets

    QVector<QString> checkbox_name;///< vector of the label texts of the QCheckBox widgets
    QVector<QCheckBox*> checkbox;///< vector of pointers to the QCheckBox widgets

    QVector<QString> double_lineedit_name;///< vector of the label texts of the double QLineEdit widgets
    QVector<QLineEdit*> double_lineedit;///< vector of pointers to the double QLineEdit widgets

    QVector<QString> string_lineedit_name;///< vector of the label texts of the string QLineEdit widgets
    QVector<QLineEdit*> string_lineedit;///< vector of pointers to the string QLineEdit widgets

    QVector<QString> filename_dialoglineedit_name;///< vector of the label texts of the filename DialogLineEdit widgets
    QVector<DialogLineEdit*> filename_dialoglineedit;///< vector of pointers to the filename DialogLineEdit widgets

    QVector<QString> directory_dialoglineedit_name;///< vector of the label texts of the directory DialogLineEdit widgets
    QVector<DialogLineEdit*> directory_dialoglineedit;///< vector of pointers to the directory DialogLineEdit widgets
  
  public:
    //constructors
    /**
     * Constructor using the (org,app) pair to determine QSettings
     * @param org organization
     * @param app application
     * @param group group
     * @param parent Parent QWidget
     */
    GuiSettingsTab(QString org, QString app, QString group, QWidget *parent = 0);

};

#endif
