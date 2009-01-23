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
