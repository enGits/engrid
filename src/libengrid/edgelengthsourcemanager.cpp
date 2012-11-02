// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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

#include "edgelengthsourcemanager.h"

#include "guiedgelengthsourcesphere.h"
#include "guiedgelengthsourcecone.h"
#include "guiedgelengthsourcebox.h"
#include "guiedgelengthsourcepipe.h"
#include "guimainwindow.h"

EdgeLengthSourceManager::EdgeLengthSourceManager()
{
  m_Sources.clear();
  m_Samples.clear();
  m_Samples.push_back(new GuiEdgeLengthSourceSphere);
  m_Samples.push_back(new GuiEdgeLengthSourceCone);
  m_Samples.push_back(new GuiEdgeLengthSourceBox);
  m_Samples.push_back(new GuiEdgeLengthSourcePipe);
  m_ListWidget = NULL;
}

EdgeLengthSourceManager::~EdgeLengthSourceManager()
{
  foreach (EdgeLengthSource* source, m_Sources) {
    //delete source;
  }
  foreach (EdgeLengthSource* sample, m_Samples) {
    //delete sample;
  }
}

void EdgeLengthSourceManager::populateListWidget()
{
  if (m_ListWidget) {
    m_ListWidget->clear();
    foreach (EdgeLengthSource *source, m_Sources) {
      m_ListWidget->addItem(source->name());
    }
  }
}

void EdgeLengthSourceManager::read()
{
  foreach (EdgeLengthSource* source, m_Sources) {
    //delete source;
  }
  m_Sources.clear();
  QString xml_text = GuiMainWindow::pointer()->getXmlSection("engrid/sources");
  QStringList lines = xml_text.split("\n");
  foreach (QString line, lines) {
    foreach (EdgeLengthSource* sample, m_Samples) {
      if (sample->read(line.trimmed())) {
        if (dynamic_cast<GuiEdgeLengthSourceSphere*>(sample)) {
          GuiEdgeLengthSourceSphere *S = new GuiEdgeLengthSourceSphere;
          S->read(line.trimmed());
          m_Sources.push_back(S);
          break;
        }
        if (dynamic_cast<GuiEdgeLengthSourceCone*>(sample)) {
          GuiEdgeLengthSourceCone *S = new GuiEdgeLengthSourceCone;
          S->read(line.trimmed());
          m_Sources.push_back(S);
          break;
        }
        if (dynamic_cast<GuiEdgeLengthSourceBox*>(sample)) {
          GuiEdgeLengthSourceBox *S = new GuiEdgeLengthSourceBox;
          S->read(line.trimmed());
          m_Sources.push_back(S);
          break;
        }
        if (dynamic_cast<GuiEdgeLengthSourcePipe*>(sample)) {
          GuiEdgeLengthSourcePipe *S = new GuiEdgeLengthSourcePipe;
          S->read(line.trimmed());
          m_Sources.push_back(S);
          break;
        }
      }
    }
  }
}

void EdgeLengthSourceManager::write()
{
  QString xml_text = "";
  foreach (EdgeLengthSource* source, m_Sources) {
    xml_text += source->write();
    xml_text += "\n";
  }
  GuiMainWindow::pointer()->setXmlSection("engrid/sources", xml_text);
}

double EdgeLengthSourceManager::minEdgeLength(vec3_t x)
{
  double L_min = 1e99;
  foreach (EdgeLengthSource* source, m_Sources) {
    double L = source->edgeLength(x);
    if (L > 0) {
      L_min = min(L, L_min);
    }
  }
  return L_min;
}

void EdgeLengthSourceManager::edit()
{
  if (m_ListWidget->currentItem()) {
    QString selected_name = m_ListWidget->currentItem()->text();
    foreach (EdgeLengthSource* source, m_Sources) {
      if (source->name() == selected_name) {
        source->config();
        break;
      }
    }
    populateListWidget();
  }
}

void EdgeLengthSourceManager::remove()
{
  if (m_ListWidget->currentItem()) {
    QList<EdgeLengthSource*> new_sources;
    QString selected_name = m_ListWidget->currentItem()->text();
    foreach (EdgeLengthSource* source, m_Sources) {
      if (source->name() == selected_name) {
        delete source;
      } else {
        new_sources.push_back(source);
      }
    }
    m_Sources = new_sources;
    populateListWidget();
  }
}

void EdgeLengthSourceManager::addSphere()
{
  QString name = "sphere" + timeStamp();
  GuiEdgeLengthSourceSphere *S = new GuiEdgeLengthSourceSphere;
  S->setName(name);
  m_Sources.append(S);
  populateListWidget();
}

void EdgeLengthSourceManager::addCone()
{
  QString name = "cone" + timeStamp();
  GuiEdgeLengthSourceCone *S = new GuiEdgeLengthSourceCone;
  S->setName(name);
  m_Sources.append(S);
  populateListWidget();
}

void EdgeLengthSourceManager::addPipe()
{
  QString name = "pipe" + timeStamp();
  GuiEdgeLengthSourcePipe *S = new GuiEdgeLengthSourcePipe;
  S->setName(name);
  m_Sources.append(S);
  populateListWidget();
}

void EdgeLengthSourceManager::addBox()
{
  QString name = "box" + timeStamp();
  GuiEdgeLengthSourceBox *S = new GuiEdgeLengthSourceBox;
  S->setName(name);
  m_Sources.append(S);
  populateListWidget();
}



