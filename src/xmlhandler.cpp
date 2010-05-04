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

#include <iostream>
using namespace std;

XmlHandler::XmlHandler(QString tagName, QObject *parent)
    : QObject(parent) {
  m_parent = parent;

  // initialise XML document
  m_TagName = tagName;
  QDomElement root = m_XmlDoc.createElement(m_TagName);
  m_XmlDoc.appendChild(root);
}

bool XmlHandler::openXml(QString file_name) {
  QFile xml_file(file_name);
  if (!xml_file.open(QIODevice::ReadOnly)) {
    QString error_message = "Failed to open xml_file " + xml_file.fileName();
    error_message += QString("\n") + "QDir::current()=" + QDir::current().absolutePath();
    error_message += QString("\n") + "QDir::currentPath()=" + QDir::currentPath();

    qWarning() << error_message;
    return(false);
  }
  if (!m_XmlDoc.setContent(&xml_file)) {
    qWarning() << tr("Error reading XML file:\n").arg(file_name) << "\n setContent failed.";
    return(false);
  }
  xml_file.close();

  // initialize m_DomNode
  resetToTopNode();

  return(true);
}

bool XmlHandler::saveXml(QString file_name) {
  QString buffer = m_XmlDoc.toString(2);
//   QString buffer = m_XmlDoc.toString(0);
  QFile xml_file(file_name);
  xml_file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream f(&xml_file);
  f << buffer << endl;
  xml_file.close();
  return(true);
}

QString XmlHandler::getXmlSection(QString name) {
  QStringList tags = name.toLower().split(tr("/"), QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  bool found = true;
  QString section_text = tr("");
  try {
    foreach(QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN(tr("error retrieving XML section '") + name + tr("'"));
      }
      if (nodes.size() == 0) {
        found = false;
        break;
      }
      if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN(tr("error retrieving XML section '") + name + tr("'"));
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

void XmlHandler::setXmlSection(QString name, QString contents) {
  QStringList tags = name.toLower().split(tr("/"), QString::SkipEmptyParts);
  QDomElement element = m_XmlDoc.documentElement();
  try {
    foreach(QString tag, tags) {
      QDomNodeList nodes = element.elementsByTagName(tag);
      if (nodes.size() > 1) {
        EG_ERR_RETURN(tr("error retrieving XML section '") + name + tr("'"));
      }
      if (nodes.size() == 0) {
        QDomElement new_element = m_XmlDoc.createElement(tag);
        element.appendChild(new_element);
        element = new_element;
      } else if (!nodes.at(0).isElement()) {
        EG_ERR_RETURN(tr("error retrieving XML section '") + name + tr("'"));
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

QStringList XmlHandler::allKeys() {
  QStringList ret;
  QDomNode dom_node;
  dom_node = m_XmlDoc.firstChild();
  parseNode(dom_node, ret, dom_node.nodeName());
  return(ret);
}

QStringList XmlHandler::childGroups() {
  QStringList ret;

  for (QDomNode sub_node = m_DomNode.firstChild(); !sub_node.isNull(); sub_node = sub_node.nextSibling()) {
    if (sub_node.nodeType() == QDomNode::ElementNode) {
      ret << sub_node.nodeName();
    }
  }

  return(ret);
}

QStringList XmlHandler::childKeys() {
  QStringList ret;
  qDebug() << "m_DomNode.nodeName()=" << m_DomNode.nodeName();
  QString output = parseNode(m_DomNode, ret, m_DomNode.nodeName());
  qDebug() << "output=" << output;
  qDebug() << "ret=" << ret;
  return(ret);
}

QString XmlHandler::group(bool absolute) {
  if (!absolute) {
    return m_DomNode.nodeName();
  } else {
    QString section;
    QString stop_node = m_XmlDoc.firstChild().nodeName();
    QDomNode parent_node = m_DomNode;
    while (!parent_node.isNull() && parent_node.nodeName() != stop_node) {
      section = parent_node.nodeName() + tr("/") + section;
      parent_node = parent_node.parentNode();
    }
    return section;
  }
}

void XmlHandler::setGroup(const QString & prefix) {
  resetToTopNode();
  beginGroup(prefix);
}

void XmlHandler::beginGroup(const QString & prefix) {
  QStringList tag_list = prefix.split(tr("/"));
  foreach(QString tag, tag_list) {
    if (!tag.isEmpty()) {
      QDomElement dom_element = m_DomNode.toElement();
      QDomNodeList dom_node_list = dom_element.elementsByTagName(tag);
      m_DomNode = dom_node_list.at(0);
    }
  }
}

void XmlHandler::endGroup() {
  if (!m_DomNode.parentNode().isNull()) m_DomNode = m_DomNode.parentNode();
}

QString XmlHandler::parseNode(const QDomNode& dom_node, QStringList& string_list, QString stop_node) {
  QString section;

  if (dom_node.nodeType() == QDomNode::TextNode) {
    QDomNode parent_node = dom_node.parentNode();
    while (!parent_node.isNull() && parent_node.nodeName() != stop_node) {
      section = parent_node.nodeName() + tr("/") + section;
      parent_node = parent_node.parentNode();
    }
    string_list << section;
  }

  for (QDomNode sub_node = dom_node.firstChild(); !sub_node.isNull(); sub_node = sub_node.nextSibling()) {
    parseNode(sub_node, string_list, stop_node);
  }

//   qDebug() << "section=" << section;
  return(section);
}

void XmlHandler::resetToTopNode() {
  m_DomNode = m_XmlDoc.firstChild();
}

void XmlHandler::resetXmlDoc() {
  m_XmlDoc.clear();
  QDomElement root = m_XmlDoc.createElement(m_TagName);
  m_XmlDoc.appendChild(root);
}
