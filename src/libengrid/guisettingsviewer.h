// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#include "guisettingstab.h"

// class QPushButton;
// class QSettings;
// class QTreeWidget;
// class QTreeWidgetItem;

// IMPORTANT:
// Setting keys can contain any Unicode characters. The Windows registry and INI files use case-insensitive keys, whereas the Carbon Preferences API on Mac OS X uses case-sensitive keys. To avoid portability problems, follow these two simple rules:
// Always refer to the same key using the same case. For example, if you refer to a key as "text fonts" in one place in your code, don't refer to it as "Text Fonts" somewhere else.
//   Avoid key names that are identical except for the case. For example, if you have a key called "MainWindow", don't try to save another key as "mainwindow".
//     Do not use slashes ('/' and '\') in key names; the backslash character is used to separate sub keys (see below). On windows '\' are converted by QSettings to '/', which makes them identical.

/**
  * Creates a QWidget with one tab per main group found in the specified QSettings file. Each of those tabs is a SettingsTab.
  */
class GuiSettingsViewer : public QDialog
{
    Q_OBJECT

  public:

    // constructors
    /**
     * Constructor using the (org,app) pair to determine QSettings
     * @param org organization
     * @param app application
     * @param parent Parent QWidget
     */
    GuiSettingsViewer(QString org, QString app, QWidget *parent = 0);

    /**
     * Constructor taking a QSettings argument to build the widget.
     * @param Set QSettings to use
     * @param parent Parent QWidget
     */
    GuiSettingsViewer(QSettings* Set, QWidget *parent = 0);

  public:
    void CreateViewer();

  private slots:
    void open();
    void save();///< save the settings
    void readSettings();

    /**  add child settings
    * @todo Delete the tabs for real. Or make sure they have the correct parent.
    */
    void addChildSettings();

  private:

    QPushButton *openButton;///< Button for the open() action
    QPushButton *saveButton;///< Button for the save() action
    QPushButton *closeButton;///< Button for the close() action

    QTabWidget tabWidget;///< a QTabWidget
    QVector<GuiSettingsTab> tabs;///< a vector containing the tabs

    QString organization;///< organization: cf QSettings documentation for more info
    QString application;///< application: cf QSettings documentation for more info
    QSettings* m_settings;///< The settings used

  public:
    /**
     * if key=value pair not found in m_settings file, write it + read key value from m_settings file and return it.
     * Version for int variables
     * @param group group
     * @param key key
     * @param value value
     */
    int getSet(QString group, QString key, int value);

    /**
     * if key=value pair not found in m_settings file, write it + read key value from m_settings file and return it.
     * Version for double variables
     * @param group group
     * @param key key
     * @param value value
     */
    double getSet(QString group, QString key, double value);

    /**
     * if key=value pair not found in m_settings file, write it + read key value from m_settings file and return it.
     * Version for bool variables
     * @param group group
     * @param key key
     * @param value value
     */
    bool getSet(QString group, QString key, bool value);

    /**
     * if key=value pair not found in m_settings file, write it + read key value from m_settings file and return it.
     * Version for path variables
     * @param group group
     * @param key key
     * @param value value
     * @param type type \n
     * type = 0 : standard string \n
     * type = 1 : filename \n
     * type = 2 : directory
     */
    QString getSet(QString group, QString key, QString value, int type);
};

#endif
