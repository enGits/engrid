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

#include "createvolumemesh.h"
#include "deletetetras.h"
#include "guimainwindow.h"
#include "updatedesiredmeshdensity.h"
#include <vtkXMLUnstructuredGridWriter.h>

CreateVolumeMesh::CreateVolumeMesh()
{
  EG_TYPENAME;
}

void CreateVolumeMesh::setTraceCells(const QVector<vtkIdType> &cells)
{
  m_TraceCells.resize(cells.size());
  qCopy(cells.begin(), cells.end(), m_TraceCells.begin());
}

void CreateVolumeMesh::getTraceCells(QVector<vtkIdType> &cells)
{
  cells.resize(m_TraceCells.size());
  qCopy(m_TraceCells.begin(), m_TraceCells.end(), cells.begin());
}

void CreateVolumeMesh::computeMeshDensity()
{
  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings").replace("\n", " ");
  if (!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaxEdgeLength;
    in >> m_MinEdgeLength;
    in >> m_GrowthFactor;
  } else {
    m_MaxEdgeLength = 1000.0;
    m_MinEdgeLength = 0.0;
    m_GrowthFactor = 1.5;
  }
  m_ELSManager.read();
  QVector<double> H(m_Grid->GetNumberOfPoints(), m_MaxEdgeLength);

  QVector<bool> fixed(m_Grid->GetNumberOfPoints(), false);
  double H_min = 1e99;
  vtkIdType id_min = -1;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    bool volume_only = true;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      if (isSurface(m_Part.n2cGG(id_node, i), m_Grid)) {
        volume_only = false;
      }
      if (!volume_only) {
        fixed[id_node] = true;
        H[id_node] = 0;
        int N = 0;
        vec3_t xi;
        m_Grid->GetPoint(id_node, xi.data());
        for (int j = 0; j < m_Part.n2nGSize(id_node); ++j) {
          if (m_Part.n2nGG(id_node, j)) {
            vec3_t xj;
            m_Grid->GetPoint(m_Part.n2nGG(id_node, j), xj.data());
            H[id_node] += (xi-xj).abs();
            ++N;
          }
        }
        if (N < 2) {
          EG_BUG;
        }
        H[id_node] /= N;
		if (H[id_node] < 0) {
	      EG_BUG;
		}
        if (H[id_node] < H_min) {
          id_min = id_node;
          H_min = H[id_node];
        }
      }
    }
  }
  if (id_min < 0) {
    EG_BUG;
  }

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    double cl_src = m_ELSManager.minEdgeLength(x);
    if (cl_src > 0) {
      if (cl_src < H[id_node]) {
        H[id_node] = cl_src;
      }
    }
  }

  QVector<bool> marked(m_Grid->GetNumberOfPoints(), false);
  marked[id_min] = true;
  bool done = false;
  while (!done) {
    done = true;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (marked[id_node] && H[id_node] <= H_min) {
        vec3_t x1;
        m_Grid->GetPoint(id_node, x1.data());
        for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node,i);
          if (!marked[id_neigh]) {
            vec3_t x2;
            m_Grid->GetPoint(id_neigh, x2.data());
            double dist = (x1 - x2).abs();
            double h = H[id_node] + (m_GrowthFactor - 1)*dist;
            if (h < 0) {
              EG_BUG;
            }
            H[id_neigh] = min(H[id_neigh], h);
            marked[id_neigh] = true;
            //H[id_neigh] += 1.0*H[id_node];
            done = false;
          }
        }
      }
    }
    H_min *= m_GrowthFactor;
  }

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    double cl_src = m_ELSManager.minEdgeLength(x);
    if (cl_src > 0) {
      if (cl_src < H[id_node]) {
        H[id_node] = cl_src;
      }
    }
  }

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x1, x2;
    m_Grid->GetPoint(id_node, x1.data());
    x2 = x1;
    for (int j = 0; j < m_Part.n2nGSize(id_node); ++j) {
      vec3_t xj;
      m_Grid->GetPoint(m_Part.n2nGG(id_node, j), xj.data());
      for (int k = 0; k < 3; ++k) {
        x1[k] = min(xj[k], x1[k]);
        x2[k] = max(xj[k], x2[k]);
      }
    }

    // Do something here!!
    EG_BUG;

  }
}


void CreateVolumeMesh::operate()
{
}

