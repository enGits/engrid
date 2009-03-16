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
#include <QTableWidgetItem>

class Cell : public QTableWidgetItem
{
public:
  Cell();
  
  QTableWidgetItem *clone() const;
  void setData(int role, const QVariant &value);
  QVariant data(int role) const;
  void setFormula(const QString &formula);
  QString formula() const;
  void setDirty();
  
private:
  QVariant value() const;
  QVariant evalExpression(const QString &str, int &pos) const;
  QVariant evalTerm(const QString &str, int &pos) const;
  QVariant evalFactor(const QString &str, int &pos) const;
  
  mutable QVariant cachedValue;
  mutable bool cacheIsDirty;
};

class SettingsSheet : public QTableWidget
{
Q_OBJECT
public:
    SettingsSheet(QWidget *parent = 0);

    ~SettingsSheet();

  bool readFile(const QString &fileName,int verbose=1);
  bool writeFile(const QString &fileName);
  void setFormula(int row, int column, const QString &formula);
  QString formula(int row, int column) const;
  Cell* cell(int row, int column) const;
    
private:
  enum { MagicNumber = 0x7F51C883 };
  
};

#endif
