// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "xmlhandler.h"

#include <QFile>
#include <QtDebug>
#include <QDir>
#include <QString>
#include <QMessageBox>
#include <QRegularExpression>

#include "engrid.h"

#include <iostream>
using namespace std;

namespace {

QString escapeXmlText(const QString &text)
{
  QString out = text;
  QRegularExpression amp_re("&(?![a-zA-Z]+;|#\\d+;|#x[0-9A-Fa-f]+;)");
  out.replace(amp_re, "&amp;");
  out.replace("<", "&lt;");
  out.replace(">", "&gt;");
  return out;
}

QString escapeXmlSections(const QString &xml, const QStringList &tags)
{
  QString out = xml;
  foreach (const QString &tag, tags) {
    const QString open = "<" + tag + ">";
    const QString close = "</" + tag + ">";
    QString rebuilt;
    int pos = 0;
    while (true) {
      int open_idx = out.indexOf(open, pos, Qt::CaseInsensitive);
      if (open_idx < 0) {
        rebuilt += out.mid(pos);
        break;
      }
      int content_start = open_idx + open.length();
      int close_idx = out.indexOf(close, content_start, Qt::CaseInsensitive);
      if (close_idx < 0) {
        rebuilt += out.mid(pos);
        break;
      }
      rebuilt += out.mid(pos, content_start - pos);
      rebuilt += escapeXmlText(out.mid(content_start, close_idx - content_start));
      rebuilt += out.mid(close_idx, close.length());
      pos = close_idx + close.length();
    }
    out = rebuilt;
  }
  return out;
}

} // namespace

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
  QString xml_text = QString::fromUtf8(xml_file.readAll());
  xml_file.close();

  QString error_msg;
  int error_line = 0;
  int error_column = 0;
  if (!m_XmlDoc.setContent(xml_text, false, &error_msg, &error_line, &error_column)) {
    QString repaired_text = escapeXmlSections(xml_text, QStringList() << "rules");
    if (repaired_text != xml_text) {
      if (!m_XmlDoc.setContent(repaired_text, false, &error_msg, &error_line, &error_column)) {
        qWarning() << tr("Error reading XML file:\n%1\n%2 (line %3, column %4)")
                      .arg(file_name)
                      .arg(error_msg)
                      .arg(error_line)
                      .arg(error_column);
        return(false);
      }
    } else {
      qWarning() << tr("Error reading XML file:\n%1\n%2 (line %3, column %4)")
                    .arg(file_name)
                    .arg(error_msg)
                    .arg(error_line)
                    .arg(error_column);
      return(false);
    }
  }

  // initialize m_DomNode
  resetToTopNode();

  return(true);
}


QString XmlHandler::nodeToString(QDomNode node, int level)
{
  QString indent_txt = "";
  for (int i = 0; i < level; ++i) {
    indent_txt += "  ";
  }
  QString txt = "";
  if (node.isElement()) {
    QDomElement elem = node.toElement();
    txt += indent_txt + "<" + elem.tagName() + ">\n";
    QDomNodeList children = node.childNodes();
    for (int i = 0; i < children.count(); ++i) {
      txt += nodeToString(children.at(i), level + 1);
    }
    txt += indent_txt + "</" + elem.tagName() + ">\n";
  }
  if (node.isText()) {
    QString raw = node.toText().data().trimmed();
    QStringList lines = raw.split('\n');
    foreach (QString line, lines) {
      txt += indent_txt + line.trimmed() + "\n";
    }
  }
  return txt;
}

bool XmlHandler::saveXml(QString file_name)
{
  QDomNodeList cases = m_XmlDoc.elementsByTagName("engridcase");
  if (cases.count() == 1) {
    QString buffer = nodeToString(cases.at(0), 0);

    QFile xml_file(file_name);
    xml_file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream f(&xml_file);
    f << buffer << endl;
    xml_file.close();
    return(true);
  }
  return false;
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
