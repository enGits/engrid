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
#ifndef SURFACEPROJECTION_H
#define SURFACEPROJECTION_H

class SurfaceProjection;

#include "egvtkobject.h"
#include "octree.h"
#include "guimainwindow.h"

class SurfaceProjection : public EgVtkObject
{

private: // data-types

  struct Triangle
  {
    vtkIdType id_a, id_b, id_c;
    vec3_t a, b, c;
    vec3_t g1, g2, g3;
    mat3_t G, GI;
    double A;
  };

  struct Edge
  {
    vec3_t a, b;
    vec3_t v;
    vec3_t na, nb;
    double La, Lb;
  };

private: // attributes

  vtkUnstructuredGrid*   m_BGrid;
  QVector<double>        m_EdgeLength;
  QVector<vtkIdType>     m_Cells;
  QVector<vtkIdType>     m_Nodes;
  QVector<vec3_t>        m_NodeNormals;
  QVector<Triangle>      m_Triangles;
  QVector<Edge>          m_Edges;
  QVector<QVector<int> > m_N2N;
  Octree                 m_OTGrid;
  QVector<double>        m_G;
  QVector<bool>          m_GSet;
  double                 m_Relax;
  double                 m_RelEdgeLength;
  double                 m_DistExp;
  double                 m_DistWeight;
  double                 m_DirWeight;
  double                 m_DirExp;
  double                 m_WeightOffset;

private: // methods

  template <class C>
  void setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells);
  void setBackgroundGrid_initOctree();
  void setBackgroundGrid_refineFromNodes();
  void setBackgroundGrid_refineFromEdges();
  void setBackgroundGrid_refineFromFaces();
  void setBackgroundGrid_initLevelSet();
  void setBackgroundGrid_computeLevelSet();

  vec3_t calcGradG(vec3_t x);
  double calcG(vec3_t x);

public: // methods

  SurfaceProjection();

  template <class C>
  void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells);
  vec3_t project(vec3_t x);
  int getNumOctreeCells() { return m_OTGrid.getNumCells(); }
  void setRelaxation(double relax) { m_Relax = relax; }
  void writeOctree(QString file_name);
  void setRelEdgeLength(double l) { m_RelEdgeLength = l; }

};

template <class C>
void SurfaceProjection::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  setBackgroundGrid_setupGrid(grid, cells);
  setBackgroundGrid_initOctree();
  setBackgroundGrid_refineFromNodes();
  setBackgroundGrid_refineFromEdges();
  setBackgroundGrid_refineFromFaces();
  setBackgroundGrid_computeLevelSet();
}

template <class C>
void SurfaceProjection::setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  QVector<vtkIdType> nodes;
  getNodesFromCells(cells, nodes, grid);
  allocateGrid(m_BGrid, cells.size(), nodes.size());
  QVector<vtkIdType> _nodes(grid->GetNumberOfPoints());
  vtkIdType id_new_node = 0;
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    grid->GetPoints()->GetPoint(id_node, x.data());
    m_BGrid->GetPoints()->SetPoint(id_new_node, x.data());
    _nodes[id_node] = id_new_node;
    ++id_new_node;
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = grid->GetCellType(id_cell);
    grid->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType new_pts[N_pts];
    for (int i = 0; i < N_pts; ++i) {
      new_pts[i] = _nodes[pts[i]];
    }
    m_BGrid->InsertNextCell(type_cell, N_pts, new_pts);
  }
  getAllCells(m_Cells, m_BGrid);
  getNodesFromCells(m_Cells, m_Nodes, m_BGrid);
  QVector<int> m_LNodes(m_Nodes.size());
  for (int i = 0; i < m_LNodes.size(); ++i) {
    m_LNodes[i] = i;
  }
  createNodeToNode(m_Cells, m_Nodes, m_LNodes, m_N2N, m_BGrid);
  m_EdgeLength.fill(1e99, m_BGrid->GetNumberOfPoints());
  foreach (vtkIdType id_node, m_Nodes) {
    vec3_t x;
    m_BGrid->GetPoints()->GetPoint(id_node, x.data());
    foreach (vtkIdType id_neigh, m_N2N[id_node]) {
      vec3_t xn;
      m_BGrid->GetPoints()->GetPoint(id_neigh, xn.data());
      m_EdgeLength[id_node] = min(m_EdgeLength[id_node], (x-xn).abs());
    }
    if (m_N2N[id_node].size() < 2) {
      EG_BUG;
    }
  }
}

#endif // SURFACEPROJECTION_H
