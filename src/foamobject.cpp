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

#include "foamobject.h"
#include <iostream>

#include "guimainwindow.h"

FoamObject::FoamObject()
{
  m_CaseDir = "";
}

int FoamObject::deleteBetween(int i, QString str1, QString str2)
{
  int i1 = m_Buffer.indexOf(str1, i);
  if (i1 != -1) {
    int i2 = m_Buffer.indexOf(str2, i1);
    if (i2 != -1) {
      m_Buffer = m_Buffer.remove(i1, i2 - i1 + str2.size());
      return i2-i1;
    }
  }
  return 0;
}

void FoamObject::stripBuffer()
{
  while (deleteBetween(0, "/*", "*/")) {};
  while (deleteBetween(0, "//", "\n")) {};
  int i = m_Buffer.indexOf("FoamFile", 0);
  if (i != -1) {
    deleteBetween(i, "{", "}");
  }
  m_Buffer.remove("FoamFile");
  m_Buffer = m_Buffer.replace("{", " ");
  m_Buffer = m_Buffer.replace("}", " ");
  m_Buffer = m_Buffer.replace("(", " ");
  m_Buffer = m_Buffer.replace(")", " ");
  m_Buffer = m_Buffer.replace(";", " ");
  m_Buffer = m_Buffer.simplified();
}

void FoamObject::readFile(QString file_name)
{
  file_name = m_CaseDir + "/" + file_name;
  QFile file(file_name);
  if (!file.open(QIODevice::ReadOnly)) {
    EG_ERR_RETURN(QString("error loading file:\n") + file_name);
  }
  QTextStream f(&file);
  m_Buffer = f.readAll();
  stripBuffer();
}

void FoamObject::buildMaps()
{
  int num_foam_nodes;
  readFile("constant/polyMesh/points");
  {
    QTextStream f(getBuffer());
    f >> num_foam_nodes;
  }
  {
    readFile("constant/polyMesh/neighbour");
    QTextStream f(getBuffer());
    int num_neigh;
    f >> num_neigh;
    int neigh = 0;
    m_FirstBoundaryFace = 0;
    while (m_FirstBoundaryFace < num_neigh) {
      f >> neigh;
      if (neigh == -1) {
        break;
      }
      ++m_FirstBoundaryFace;
    }
  }

  m_VolToSurfMap.fill(-1, num_foam_nodes);
  readFile("constant/polyMesh/faces");
  int num_surf_nodes = 0;
  int max_node = 0;
  {
    int num_foam_faces;
    QTextStream f(getBuffer());
    f >> num_foam_faces;
    for (int i = 0; i < num_foam_faces; ++i) {
      int num_nodes;
      f >> num_nodes;
      for (int j = 0; j < num_nodes; ++j) {
        int node;
        f >> node;
        max_node = max(node, max_node);
        if (i >= m_FirstBoundaryFace) {
          if (m_VolToSurfMap[node] == -1) {
            m_VolToSurfMap[node] = num_surf_nodes;
            ++num_surf_nodes;
          }
        }
      }
    }
  }
  m_SurfToVolMap.fill(-1, num_surf_nodes);
  for (int i = 0; i < m_VolToSurfMap.size(); ++i) {
    if (m_VolToSurfMap[i] != -1) {
      if (m_VolToSurfMap[i] > m_SurfToVolMap.size()) {
        EG_BUG;
      }
      m_SurfToVolMap[m_VolToSurfMap[i]] = i;
    }
  }
  for (int i = 0; i < m_SurfToVolMap.size(); ++i) {
    if (m_SurfToVolMap[i] == -1) {
      EG_BUG;
    }
  }
}

void FoamObject::setCaseDir(QString case_dir)
{
  m_CaseDir = case_dir;
  GuiMainWindow::pointer()->setXmlSection("openfoam/CaseDir",m_CaseDir);
}
