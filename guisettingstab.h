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
  QVector<QString> spinbox_name;
  QVector<QSpinBox*> spinbox;
  
  QVector<QString> checkbox_name;
  QVector<QCheckBox*> checkbox;
  
  QVector<QString> lineedit_name;
  QVector<QLineEdit*> lineedit;
  
public:
	//constructors
    SettingsTab(QWidget *parent = 0);
    SettingsTab(QString org,QString app,QString group,QWidget *parent = 0);
};

#endif
