//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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

#include <vtkDelaunay2D.h>
#include <vtkGeometryFilter.h>
#include <vtkCellArray.h>

#include "fillplane.h"
#include "vtkEgPolyDataToUnstructuredGridFilter.h"
#include "guimainwindow.h"
#include "correctsurfaceorientation.h"

FillPlane::FillPlane()
{
  m_InverseDirection = false;
  m_BC = 0;
}

vec3_t FillPlane::toPlane(vec3_t x)
{
  x -= m_X0;
  x -= (x*m_N)*m_N;
  return vec3_t(x*m_G1, x*m_G2, 0);
}

vec3_t FillPlane::fromPlane(vec3_t x)
{
  return m_X0 + x[0]*m_G1 + x[1]*m_G2;
}

bool FillPlane::isWithinTolerance(vec3_t x)
{
  if (fabs((x - m_X0)*m_N) < m_Tol) {
    return true;
  }
  return false;
}

void FillPlane::createEdgesOnPlane(vtkUnstructuredGrid *edge_grid)
{
  vtkIdType num_edges = 0;
  vtkIdType num_nodes = 0;
  QVector<vtkIdType> node_map(m_Grid->GetNumberOfPoints(), -1);
  for (vtkIdType id_node1 = 0; id_node1 < m_Grid->GetNumberOfPoints(); ++id_node1) {
    vec3_t x;
    m_Grid->GetPoint(id_node1, x.data());
    if (isWithinTolerance(x)) {
      node_map[id_node1] = num_nodes;
      ++num_nodes;
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        if (node_map[id_node2] != -1) {
          ++num_edges;
        }
      }
    }
  }
  m_NodeMap.resize(num_nodes);
  allocateGrid(edge_grid, 2*num_edges, num_nodes, false);
  for (vtkIdType id_node1 = 0; id_node1 < m_Grid->GetNumberOfPoints(); ++id_node1) {
    if (node_map[id_node1] != -1) {
      m_NodeMap[node_map[id_node1]] = id_node1;
      vec3_t x;
      m_Grid->GetPoint(id_node1, x.data());
      edge_grid->GetPoints()->SetPoint(node_map[id_node1], x.data());
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        if (node_map[id_node2] != -1 && id_node2 > id_node1) {
          vtkIdType pts[2];
          pts[0] = node_map[id_node1];
          pts[1] = node_map[id_node2];
          edge_grid->InsertNextCell(VTK_LINE, 2, pts);
        }
      }
    }
  }
}

void FillPlane::closeLoops(vtkUnstructuredGrid *edge_grid)
{
  bool done = false;
  while (!done) {
    QList<vtkIdType> end_nodes;
    QVector<int> count(edge_grid->GetNumberOfPoints(), 0);
    for (vtkIdType id_edge = 0; id_edge < edge_grid->GetNumberOfCells(); ++id_edge) {
      vtkIdType num_pts, *pts;
      edge_grid->GetCellPoints(id_edge, num_pts, pts);
      for (int i = 0; i < num_pts; ++i) {
        ++count[pts[i]];
      }
    }
    for (vtkIdType id_node = 0; id_node < edge_grid->GetNumberOfPoints(); ++id_node) {
      if (count[id_node] == 0) {
        EG_ERR_RETURN("unable to fill plane(s)");
      }
      if (count[id_node] > 2) {
        EG_ERR_RETURN("unable to fill plane(s)");
      }
      if (count[id_node] == 1) {
        end_nodes.append(id_node);
      }
    }
    if (end_nodes.size() % 2 != 0) {
      EG_ERR_RETURN("unable to fill plane(s)");
    }
    if (end_nodes.size() > 0) {
      double dist_min = EG_LARGE_REAL;
      vtkIdType id_fill1 = -1;
      vtkIdType id_fill2 = -1;
      foreach (vtkIdType id_node1, end_nodes) {
        vec3_t n1 = m_Part.globalNormal(m_NodeMap[id_node1]);
        vec3_t x1;
        edge_grid->GetPoint(id_node1, x1.data());
        foreach (vtkIdType id_node2, end_nodes) {
          if (id_node1 != id_node2) {
            vec3_t n2 = m_Part.globalNormal(m_NodeMap[id_node2]);
            vec3_t x2;
            edge_grid->GetPoint(id_node2, x2.data());
            if (n1*n2 < -0.5) {
              double dist = (x1 - x2).abs();
              if (dist < dist_min) {
                dist_min = dist;
                id_fill1 = id_node1;
                id_fill2 = id_node2;
              }
            }
          }
        }
      }
      if (id_fill1 == -1 || id_fill2 == -1) {
        EG_ERR_RETURN("unable to fill plane(s)");
      }
      vtkIdType pts[2];
      pts[0] = id_fill1;
      pts[1] = id_fill2;
      edge_grid->InsertNextCell(VTK_LINE, 2, pts);
    } else {
      done = true;
    }
  }
}

void FillPlane::gridToPlane(vtkUnstructuredGrid *edge_grid)
{
  for (vtkIdType id_node = 0; id_node < edge_grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    edge_grid->GetPoint(id_node, x.data());
    x = toPlane(x);
    edge_grid->GetPoints()->SetPoint(id_node, x.data());
  }
}

void FillPlane::gridFromPlane(vtkUnstructuredGrid *edge_grid)
{
  for (vtkIdType id_node = 0; id_node < edge_grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    edge_grid->GetPoint(id_node, x.data());
    x = fromPlane(x);
    edge_grid->GetPoints()->SetPoint(id_node, x.data());
  }
}

void FillPlane::triangulate(vtkPolyData *edge_pdata, vtkUnstructuredGrid *tri_grid)
{
  EG_VTKSP(vtkDelaunay2D, delaunay);
  EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, pdata2grid);
  delaunay->SetInput(edge_pdata);
  delaunay->SetSource(edge_pdata);
  delaunay->Update();
  pdata2grid->SetInput(delaunay->GetOutput());
  pdata2grid->Update();
  makeCopy(pdata2grid->GetOutput(), tri_grid);
  createBasicFields(tri_grid, tri_grid->GetNumberOfCells(), tri_grid->GetNumberOfPoints());
  QSet<int> bcs = getAllBoundaryCodes(m_Grid);
  m_BC = 0;
  foreach (int bc, bcs) {
    m_BC = max(m_BC, bc);
  }
  ++m_BC;
  EG_VTKDCC(vtkIntArray, cell_code, tri_grid, "cell_code");
  for (vtkIdType id_face = 0; id_face < tri_grid->GetNumberOfCells(); ++id_face) {
    cell_code->SetValue(id_face, m_BC);
  }
}

void FillPlane::order(vtkUnstructuredGrid *edge_grid, vtkPolyData *edge_pdata)
{
  QVector<QVector<vtkIdType> > edges(edge_grid->GetNumberOfCells(), QVector<vtkIdType>(2));
  for (vtkIdType id_edge = 0; id_edge < edge_grid->GetNumberOfCells(); ++id_edge) {
    vtkIdType num_pts, *pts;
    edge_grid->GetCellPoints(id_edge, num_pts, pts);
    if (num_pts != 2) {
      EG_BUG;
    }
    edges[id_edge][0] = pts[0];
    edges[id_edge][1] = pts[1];
  }
  for (int i = 1; i < edges.size(); ++i) {
    for (int j = i; j < edges.size(); ++j) {
      QVector<vtkIdType> tmp_edge = edges[j];
      edges[j] = edges[i];
      if (tmp_edge[0] == edges[i-1][1]) {
        edges[i][0] = tmp_edge[0];
        edges[i][1] = tmp_edge[1];
        break;
      }
      if (tmp_edge[1] == edges[i-1][1]) {
        edges[i][0] = tmp_edge[1];
        edges[i][1] = tmp_edge[0];
        break;
      }
      edges[j] = tmp_edge;
    }
  }
  QList<vtkIdType> poly_nodes;
  for (vtkIdType id_edge = 0; id_edge < edge_grid->GetNumberOfCells(); ++id_edge) {
    poly_nodes.append(edges[id_edge][0]);
  }
  orderGeometrically(edge_grid, poly_nodes);
  EG_VTKSP(vtkPoints, points);
  EG_VTKSP(vtkCellArray, polys);
  foreach (vtkIdType id_node, poly_nodes) {
    vec3_t x;
    edge_grid->GetPoint(id_node, x.data());
    points->InsertNextPoint(x.data());
  }
  EG_VTKSP(vtkIdList, pts);
  pts->SetNumberOfIds(poly_nodes.size());
  for (vtkIdType i = 0; i < pts->GetNumberOfIds(); ++i) {
    pts->SetId(i, i);
  }
  polys->InsertNextCell(pts);
  edge_pdata->SetPoints(points);
  edge_pdata->SetPolys(polys);
}

void FillPlane::orderGeometrically(vtkUnstructuredGrid* edge_grid, QList<vtkIdType>& poly_nodes)
{
  if (poly_nodes.size() < 3) {
    return;
  }
  vec3_t x_centre(0,0,0);
  int N = 0;
  foreach (vtkIdType id_node, poly_nodes) {
    vec3_t x;
    edge_grid->GetPoint(id_node, x.data());
    x_centre += x;
    ++N;
  }
  x_centre *= 1.0/N;
  double scale_sum = 0;
  {
    bool first = true;
    vec3_t x1;
    foreach (vtkIdType id_node, poly_nodes) {
      vec3_t x;
      edge_grid->GetPoint(id_node, x.data());
      if (first) {
        x1 = x;
        first = false;
      } else {
        vec3_t x2 = x;
        vec3_t u = x2 - x1;
        vec3_t v = GeometryTools::rotate(u, m_N, deg2rad(90));
        vec3_t c = x_centre - 0.5*(x1 + x2);
        c.normalise();
        scale_sum += c*v;
        x1 = x;
      }
    }
  }
  bool invert = false;
  if (scale_sum < 0) {
    if (!m_InverseDirection) {
      invert = true;
    }
  } else {
    if (m_InverseDirection) {
      invert = true;
    }
  }
  if (invert) {
    QList<vtkIdType> reverted;
    foreach (vtkIdType id_node, poly_nodes) {
      reverted.prepend(id_node);
    }
    poly_nodes = reverted;
  }
}

void FillPlane::operate()
{
  m_G1 = GeometryTools::orthogonalVector(m_N);
  m_G2 = m_N.cross(m_G1);
  m_N.normalise();
  m_G1.normalise();
  m_G2.normalise();
  EG_VTKSP(vtkUnstructuredGrid, edge_grid);
  EG_VTKSP(vtkUnstructuredGrid, tri_grid);
  EG_VTKSP(vtkPolyData, edge_pdata);
  createEdgesOnPlane(edge_grid);
  gridToPlane(edge_grid);
  closeLoops(edge_grid);
  order(edge_grid, edge_pdata);
  triangulate(edge_pdata, tri_grid);
  gridFromPlane(tri_grid);
  writeGrid(tri_grid, "filltest");
  MeshPartition tri_part(tri_grid, true);
  m_Part.addPartition(tri_part, m_Tol);
  m_Grid->Modified();
  CorrectSurfaceOrientation corr_surf;
  corr_surf();
  GuiMainWindow::pointer()->updateBoundaryCodes(true);
}


