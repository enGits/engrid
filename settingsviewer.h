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
  
private:
    void readSettings();
    void addChildSettings();

    QTreeWidget *treeWidget;
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
