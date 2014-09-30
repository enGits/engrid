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
// 
#include "geometrysmoother.h"
#include "geometrytools.h"

GeometrySmoother::GeometrySmoother()
{
  m_CornerAngle      = GeometryTools::deg2rad(20.0);
  m_RelaxationFactor = 0.1;
  m_NumIterations    = 10;
}

void GeometrySmoother::operate()
{
  // copy all node coordinates into a QVector
  QVector<vec3_t> x(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_Grid->GetPoint(id_node, x[id_node].data());
  }

  // compute "snap" points
  int N = 0;
  QVector<QList<vtkIdType> > snap_points(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    QList<vtkIdType> boundary_neighbours;
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vtkIdType id_neigh = m_Part.n2nGG(id_node, i);
      QList<vtkIdType> edge_faces;
      m_Part.getEdgeFaces(id_node, id_neigh, edge_faces);
      if (edge_faces.size() > 2) {
        EG_BUG;
      }
      else if (edge_faces.size() < 1) {
        EG_BUG;
      }
      else if (edge_faces.size() == 1) {
        boundary_neighbours << id_neigh;
      }
      else {
        snap_points[id_node] << id_neigh;
      }
    }
    if (boundary_neighbours.size() > 2) {
      EG_BUG;
    }
    if (boundary_neighbours.size() == 2) {      
      vec3_t u = x[id_node] - x[boundary_neighbours[0]];
      vec3_t v = x[boundary_neighbours[1]] - x[id_node];
      vec3_t n = m_Part.globalNormal(id_node);
      u -= (n*u)*n;
      v -= (n*v)*n;
      if (GeometryTools::angle(u, v) > m_CornerAngle) {
        snap_points[id_node].clear();
        ++N;
      }
      else {
        snap_points[id_node] = boundary_neighbours;
      }
    }
  }

  // the smoothing process starts here
  for (int iter = 0; iter < m_NumIterations; ++iter) {
    QVector<vec3_t> x_new = x;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (snap_points[id_node].size()) {
        double w = m_RelaxationFactor/snap_points[id_node].size();
        if (snap_points.size() == 2) {
          w *= 0.3333;
        }
        foreach (vtkIdType id_snap, snap_points[id_node]) {
          x_new[id_node] += w*(x[id_snap] - x[id_node]);
        }
      }
    }
    x = x_new;
  }

  // copy node coordinates back to the VTK grid
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_Grid->GetPoints()->SetPoint(id_node, x[id_node].data());
  }
}
