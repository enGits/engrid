#ifndef SETTINGSTAB_H
#define SETTINGSTAB_H

#include <QWidget>
#include <QString>
#include <QtGui>
#include <QVector>

/**
  * Creates a QWidget listing all key/value pairs contained in the group "group" of the QSettings file corresponding to the (org,app) pair.
  * integers appear in spinboxes
  * doubles appear in line edit boxes
  * booleans appear in checkboxes
  */
class SettingsTab : public QWidget
{
Q_OBJECT
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
    SettingsTab(QString org,QString app,QString group,QWidget *parent = 0);
};

#endif
