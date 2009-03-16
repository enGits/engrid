//
// C++ Implementation: settingssheet
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <QFile>
#include <QMessageBox>
#include <QApplication>
#include <QtGui>
#include "settingssheet.h"

SettingsSheet::SettingsSheet(QWidget *parent)
 : QTableWidget(parent)
{
}


SettingsSheet::~SettingsSheet()
{
}

Cell *SettingsSheet::cell(int row, int column) const
{
  return static_cast<Cell *>(item(row, column));
}

bool SettingsSheet::readFile(const QString &fileName)
{
  QFile file(fileName);
  if (!file.open(QIODevice::ReadOnly)) {
    QMessageBox::warning(this, tr("SettingsSheet"),
                         tr("Cannot read file %1:\n%2.")
                         .arg(file.fileName())
                         .arg(file.errorString()));
    return false;
  }
  
  QDataStream in(&file);
  in.setVersion(QDataStream::Qt_4_1);
  
  quint32 magic;
  in >> magic;
  if (magic != MagicNumber) {
    QMessageBox::warning(this, tr("SettingsSheet"),
                         tr("The file is not a SettingsSheet file."));
    return false;
  }
  
  clear();
  
  quint16 row;
  quint16 column;
  QString str;
  
  QApplication::setOverrideCursor(Qt::WaitCursor);
  while (!in.atEnd()) {
    in >> row >> column >> str;
    setFormula(row, column, str);
  }
  QApplication::restoreOverrideCursor();
  return true;
}

bool SettingsSheet::writeFile(const QString &fileName)
{
  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly)) {
    QMessageBox::warning(this, tr("SettingsSheet"),
                         tr("Cannot write file %1:\n%2.")
                         .arg(file.fileName())
                         .arg(file.errorString()));
    return false;
  }
  
  QDataStream out(&file);
  out.setVersion(QDataStream::Qt_4_1);
  
  out << quint32(MagicNumber);
  
  QApplication::setOverrideCursor(Qt::WaitCursor);
  int RowCount=this->rowCount();
  int ColumnCount=this->columnCount();
  for (int row = 0; row < RowCount; ++row) {
    for (int column = 0; column < ColumnCount; ++column) {
      QString str = formula(row, column);
      if (!str.isEmpty())
        out << quint16(row) << quint16(column) << str;
    }
  }
  QApplication::restoreOverrideCursor();
  return true;
}

void SettingsSheet::setFormula(int row, int column,
                             const QString &formula)
{
  Cell *c = cell(row, column);
  if (!c) {
    c = new Cell;
    setItem(row, column, c);
  }
  c->setFormula(formula);
}

QString SettingsSheet::formula(int row, int column) const
{
  Cell *c = cell(row, column);
  if (c) {
    return c->formula();
  } else {
    return "";
  }
}

Cell::Cell()
{
  setDirty();
}

QTableWidgetItem *Cell::clone() const
{
  return new Cell(*this);
}

void Cell::setData(int role, const QVariant &value)
{
  QTableWidgetItem::setData(role, value);
  if (role == Qt::EditRole)
    setDirty();
}

QVariant Cell::data(int role) const
{
  if (role == Qt::DisplayRole) {
    if (value().isValid()) {
      return value().toString();
    } else {
      return "####";
    }
  } else if (role == Qt::TextAlignmentRole) {
    if (value().type() == QVariant::String) {
      return int(Qt::AlignLeft | Qt::AlignVCenter);
    } else {
      return int(Qt::AlignRight | Qt::AlignVCenter);
    }
  } else {
    return QTableWidgetItem::data(role);
  }
}

void Cell::setFormula(const QString &formula)
{
  setData(Qt::EditRole, formula);
}

QString Cell::formula() const
{
  return data(Qt::EditRole).toString();
}

void Cell::setDirty()
{
  cacheIsDirty = true;
}

const QVariant Invalid;

QVariant Cell::value() const
{
  if (cacheIsDirty) {
    cacheIsDirty = false;
    
    QString formulaStr = formula();
    if (formulaStr.startsWith('\'')) {
      cachedValue = formulaStr.mid(1);
    } else if (formulaStr.startsWith('=')) {
      cachedValue = Invalid;
      QString expr = formulaStr.mid(1);
      expr.replace(" ", "");
      expr.append(QChar::Null);
      
      int pos = 0;
      cachedValue = evalExpression(expr, pos);
      if (expr[pos] != QChar::Null)
        cachedValue = Invalid;
    } else {
      bool ok;
      double d = formulaStr.toDouble(&ok);
      if (ok) {
        cachedValue = d;
      } else {
        cachedValue = formulaStr;
      }
    }
  }
  return cachedValue;
}

QVariant Cell::evalExpression(const QString &str, int &pos) const
{
  QVariant result = evalTerm(str, pos);
  while (str[pos] != QChar::Null) {
    QChar op = str[pos];
    if (op != '+' && op != '-')
      return result;
    ++pos;
    
    QVariant term = evalTerm(str, pos);
    if (result.type() == QVariant::Double
        && term.type() == QVariant::Double) {
          if (op == '+') {
            result = result.toDouble() + term.toDouble();
          } else {
            result = result.toDouble() - term.toDouble();
          }
        } else {
          result = Invalid;
        }
  }
  return result;
}

QVariant Cell::evalTerm(const QString &str, int &pos) const
{
  QVariant result = evalFactor(str, pos);
  while (str[pos] != QChar::Null) {
    QChar op = str[pos];
    if (op != '*' && op != '/')
      return result;
    ++pos;
    
    QVariant factor = evalFactor(str, pos);
    if (result.type() == QVariant::Double
        && factor.type() == QVariant::Double) {
          if (op == '*') {
            result = result.toDouble() * factor.toDouble();
          } else {
            if (factor.toDouble() == 0.0) {
              result = Invalid;
            } else {
              result = result.toDouble() / factor.toDouble();
            }
          }
        } else {
          result = Invalid;
        }
  }
  return result;
}

QVariant Cell::evalFactor(const QString &str, int &pos) const
{
  QVariant result;
  bool negative = false;
  
  if (str[pos] == '-') {
    negative = true;
    ++pos;
  }
  
  if (str[pos] == '(') {
    ++pos;
    result = evalExpression(str, pos);
    if (str[pos] != ')')
      result = Invalid;
    ++pos;
  } else {
    QRegExp regExp("[A-Za-z][1-9][0-9]{0,2}");
    QString token;
    
    while (str[pos].isLetterOrNumber() || str[pos] == '.') {
      token += str[pos];
      ++pos;
    }
    
    if (regExp.exactMatch(token)) {
      int column = token[0].toUpper().unicode() - 'A';
      int row = token.mid(1).toInt() - 1;
      
      Cell *c = static_cast<Cell *>(
                                     tableWidget()->item(row, column));
      if (c) {
        result = c->value();
      } else {
        result = 0.0;
      }
    } else {
      bool ok;
      result = token.toDouble(&ok);
      if (!ok)
        result = Invalid;
    }
  }
  
  if (negative) {
    if (result.type() == QVariant::Double) {
      result = -result.toDouble();
    } else {
      result = Invalid;
    }
  }
  return result;
}
