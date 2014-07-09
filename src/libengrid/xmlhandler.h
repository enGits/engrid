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
#ifndef XMLHANDLER_H
#define XMLHANDLER_H

#include <QDomDocument>
#include <QObject>
#include <QMessageBox>
#include <QtDebug>

/** A class to simplify the handling of XML files. */
class XmlHandler : public QObject {
    Q_OBJECT
  private:
    QDomDocument m_XmlDoc; ///< XML document describing the templates to use
    QObject *m_parent;///< Parent widget (for message boxes)
    QDomNode m_DomNode;///< current node in the DOM tree
    QString m_TagName;

  public:
    QDomNode getDomNode() { return m_DomNode; }///< Returns m_DomNode

  public:
    XmlHandler(QString tagName, QObject *parent = 0);///< Constructor
    bool openXml(QString file_name);///< Open XML file
    bool saveXml(QString file_name);///< Save XML file
    QString getXmlSection(QString name);///< get contents of XML section
    void setXmlSection(QString name, QString contents);///< set contents of XML section
    void resetXmlDoc();///< Initialize or reset m_XmlDoc
    QString getBuffer(int indent = 1) { return m_XmlDoc.toString(indent); }
    QDomDocument* getXmlDoc() {return &m_XmlDoc;}

  public:
    /// Returns a list of all keys, including subkeys, that can be read using the XmlHandler object.
    QStringList allKeys();
    /// Returns a list of all key top-level groups that contain keys that can be read using the XmlHandler object.
    QStringList childGroups();
    /// Returns a list of all top-level keys that can be read using the XmlHandler object.
    QStringList childKeys();
    /// Returns the current group.
    QString group(bool absolute = true);
    ///Appends prefix to the current group.
    void beginGroup(const QString & prefix);
    /// Resets the group to what it was before the corresponding beginGroup() call.
    void endGroup();

  public:
    void resetToTopNode();///< Reset to top node
    void setGroup(const QString & prefix);///< Set current group

  public:
    /** Parses a QDomNode for text nodes.
     * @param dom_node QDomNode to parse
     * @param string_list QStringList in which to store the paths containing text
     * @param stop_node node at which to stop working up the path. ex: if the full path of a located text node is "#document/A/B/C/D/" and stop_node="B" then "C/D/" will be put into string_list (and eventually returned if dom_node is a text node)
     * @return If dom_node is a text node, returns the path to it, otherwise returns an empty string.
     */
    QString parseNode(const QDomNode& dom_node, QStringList& string_list, QString stop_node);
};

#endif
