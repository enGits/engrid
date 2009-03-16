//
// C++ Interface: settingssheet
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SETTINGSSHEET_H
#define SETTINGSSHEET_H

#include <QTableWidget>

/**
	@author Mike Taverne <mtaverne@engits.com>
*/
class SettingsSheet : public QTableWidget
{
Q_OBJECT
public:
    SettingsSheet(QWidget *parent = 0);

    ~SettingsSheet();

};

#endif
