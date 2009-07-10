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
#include "xmlhandler.h"

#include <QFile>
#include <QtDebug>
#include <QDir>
#include <QString>
#include <QMessageBox>

#include "engrid.h"

XmlHandler::XmlHandler(QObject *parent)
: QObject(parent)
{
  m_parent = parent;
}


XmlHandler::~XmlHandler()
{
}

void XmlHandler::openXml(QString file_name)
{
  qDebug()<<"Opening "<<file_name;
  QFile xml_file(file_name);
  if (!xml_file.open(QIODevice::ReadOnly)) {
    qWarning()<<"Failed to open xml_file "<<xml_file.fileName();
    qWarning()<<"QDir::current()="<<QDir::current();
    qWarning()<<"QDir::currentPath()="<<QDir::currentPath();
    EG_BUG;
  }
  if (!m_XmlDoc.setContent(&xml_file)) {
    QMessageBox::critical((QWidget*)m_parent, "Open failed", "Error reading enGrid case file:\n" + file_name);
  }
  xml_file.close();
}

void XmlHandler::saveXml(QString file_name)
{
  QString buffer = m_XmlDoc.toString(2);
  QFile xml_file(file_name);
  xml_file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&xml_file);
  f << buffer << endl;
}

QString XmlHandler::getXmlSection(QString name)
{
  QStringList tags = name.toLower().split("/", QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  bool found = true;
  QString section_text = "";
  try {
    foreach (QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      if (nodes.size() == 0) {
        found = false;
        break;
      }
      if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      element = nodes.at(0).toElement();
    }
  } catch (Error err) {
    err.display();
  }
  if (found) {
    section_text = element.text();
  }
  return section_text;
}

void XmlHandler::setXmlSection(QString name, QString contents)
{
  QStringList tags = name.toLower().split("/", QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  try {
    foreach (QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      }
      if (nodes.size() == 0) {
        QDomElement new_element = m_XmlDoc.createElement(tag);
        element.appendChild(new_element);
        element = new_element;
      } else if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN("error retrieving XML section '" + name + "'");
      } else {
        element = nodes.at(0).toElement();
      }
    }
    while (element.hasChildNodes()) {
      element.removeChild(element.firstChild());
    }
    QDomText text_node = m_XmlDoc.createTextNode(contents);
    element.appendChild(text_node);
  } catch (Error err) {
    err.display();
  }
}
