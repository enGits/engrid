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
#ifndef XMLHANDLER_H
#define XMLHANDLER_H

#include <QDomDocument>
#include <QObject>
#include <QMessageBox>
#include <QtDebug>

class XmlHandler : public QObject {
  Q_OBJECT
private:
  QDomDocument m_XmlDoc; ///< XML document describing the templates to use
  QObject *m_parent;
public:
  XmlHandler(QObject *parent = 0);
  ~XmlHandler();
  void openXml(QString file_name);
  void saveXml(QString file_name);
  QString getXmlSection(QString name);
  void setXmlSection(QString name, QString contents);
};

#endif
