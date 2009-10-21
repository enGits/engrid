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

  bool readFile(int verbose = 1);
  void setFormula(int row, int column, const QString &formula);
  QString formula(int row, int column) const;
  Cell* cell(int row, int column) const;

public slots:

  void writeFile();

};

#endif
