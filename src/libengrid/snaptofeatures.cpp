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

#include "snaptofeatures.h"
#include "guimainwindow.h"

SnapToFeatures::SnapToFeatures()
{
}

void SnapToFeatures::snapIteration()
{
  EG_VTKDCN(vtkCharArray_t, node_type, m_Grid, "node_type");
  QVector<bool> node_blocked = m_NodeSnapped;
  foreach (vtkIdType id_node, m_Part.getNodes()) {
    if (!node_blocked[id_node]) {
      --m_NodesToSnap;
      if (node_type->GetValue(id_node) == EG_SIMPLE_VERTEX && m_Part.n2bcGSize(id_node) == 1) {
        CadInterface *cad = GuiMainWindow::pointer()->getCadInterface(m_Part.n2bcG(id_node, 0));
        vec3_t x;
        m_Grid->GetPoint(id_node, x.data());
        vec3_t x_snap = cad->snapToEdge(x);
        QVector<vec3_t> normals(m_Part.n2cGSize(id_node));
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          normals[i] = cellNormal(m_Grid, m_Part.n2cGG(id_node, i));
        }
        m_Grid->GetPoints()->SetPoint(id_node, x_snap.data());
        bool no_feature = false;
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          vec3_t n = cellNormal(m_Grid, m_Part.n2cGG(id_node, i));
          if (GeometryTools::angle(normals[i], n) > 0.5*M_PI) {
            no_feature = true;
            break;
          }
        }
        node_blocked[id_node] = true;
        m_NodeSnapped[id_node] = true;
        if (no_feature) {
          m_Grid->GetPoints()->SetPoint(id_node, x.data());
        } else {
          for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
            vtkIdType id_neigh = m_Part.n2nGSize(id_node);
            node_blocked[id_neigh] = true;
            m_NodeSnapped[id_neigh] = true;
            --m_NodesToSnap;
          }
        }
      }
    }
  }
}

void SnapToFeatures::operate()
{
  m_NodeSnapped.fill(false, m_Grid->GetNumberOfPoints());
  m_NodesToSnap = m_Part.getNumberOfNodes();
  while (m_NodesToSnap > 0) {
    snapIteration();
  }
}
