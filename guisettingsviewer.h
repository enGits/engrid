#ifndef SETTINGSVIEWER_H
#define SETTINGSVIEWER_H

#include <QDialog>
#include <QtGui>

#include "guisettingstab.h"

// class QPushButton;
// class QSettings;
// class QTreeWidget;
// class QTreeWidgetItem;

/**
  * Creates a QWidget with one tab per main group found in the specified QSettingsm file. each of those tabs is a SettingsTab .
  */
class SettingsViewer : public QDialog
{
    Q_OBJECT

public:
	//constructors
  /**
   * Constructor using the (org,app) pair to determine QSettings
   * @param org organization
   * @param app application
   * @param group group
   * @param parent Parent QWidget
   */
  SettingsViewer(QString org, QString app,QWidget *parent = 0);
  /**
   * Constructor taking a QSettings argument to build the widget.
   * @param Set QSettings to use
   * @param group group
   * @param parent Parent QWidget
   */
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
