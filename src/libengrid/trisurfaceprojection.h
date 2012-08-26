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
#ifndef TRISURFACEPROJECTION_H
#define TRISURFACEPROJECTION_H

class TriSurfaceProjection;

#include "egvtkobject.h"
#include "guimainwindow.h"
#include "geometrytools.h"
#include "vtkCharArray.h"
#include "surfaceoperation.h"
#include "surfaceprojection.h"
#include "triangle.h"
#include "facefinder.h"

class TriSurfaceProjection : public SurfaceProjection
{

private: // data-types

  struct Edge
  {
    vec3_t a, b;
    vec3_t v;
    vec3_t na, nb;
    double La, Lb;
  };

private: // methods

  template <class C>
  void setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells); ///< copy the cells from grid to m_BGrid

protected: // attributes

  vtkUnstructuredGrid*      m_BGrid; ///< the background grid defining the geometry
  MeshPartition             m_BPart;
  QVector<double>           m_EdgeLength;
  QVector<vtkIdType>        m_Cells;
  QVector<vtkIdType>        m_Nodes;
  QVector<vec3_t>           m_NodeNormals; ///< The surface normal at each node of m_BGrid
  QVector<Triangle>         m_Triangles; ///< All triangles of m_BGrid. One for each triangle cell of m_BGrid.
  QVector<double>           m_Radius; ///< Surface radius for mesh resolution.
  QVector<QVector<int> >    m_N2N;
  double                    m_CritDistance;
  QMap<vtkIdType,vtkIdType> m_Pindex;
  FaceFinder                m_FaceFinder;
  bool                      m_RestrictToTriangle;
  vtkIdType                 m_LastProjTriangle;

protected: // static attributes

  static vtkIdType m_LastPindex;

protected: // methods

  virtual void   updateBackgroundGridInfo();      ///< Set up the background grid (triangles, bezier triangles, etc)

  void      searchNewTriangle(vec3_t xp, vtkIdType &id_tri, vec3_t &x_proj, vec3_t &r_proj, bool neigh_mode, bool &on_triangle);
  vtkIdType getProjTriangle(vtkIdType id_node);
  void      setProjTriangle(vtkIdType id_node, vtkIdType proj_triangle);
  void      computeSurfaceCurvature();

public: // methods

  TriSurfaceProjection();
  ~TriSurfaceProjection();
  
  template <class C> void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells); ///< Set the background grid to use + set it up

  virtual vec3_t    project(vec3_t x, vtkIdType id_node = -1, bool correct_curvature = false, vec3_t v = vec3_t(0,0,0));
  virtual double    getRadius(vtkIdType id_node);
  virtual vec3_t    correctCurvature(vtkIdType proj_triangle, vec3_t x);
  virtual vec3_t    lastProjNormal() { return GeometryTools::cellNormal(m_BGrid, m_LastProjTriangle); }
  virtual vtkIdType lastProjTriangle() { return m_LastProjTriangle; }
  virtual bool      lastProjFailed() { return false; }


  vtkUnstructuredGrid* getBGrid() { return m_BGrid; }

public: // static methods

  static void resetPindex() { m_LastPindex = 0; }

};


template <class C>
void TriSurfaceProjection::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  setBackgroundGrid_setupGrid(grid, cells);
  updateBackgroundGridInfo();
  setForegroundGrid(grid);
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_BGrid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      setProjTriangle(pts[i], id_cell);
    }
  }
  m_FaceFinder.setMaxNumFaces(10);
  m_FaceFinder.setGrid(m_BGrid);
}

template <class C>
void TriSurfaceProjection::setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  QVector<vtkIdType> nodes;
  getNodesFromCells(cells, nodes, grid);
  int num_new_cells = cells.size();
  foreach (vtkIdType id_cell, cells) {
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_QUAD) {
      ++num_new_cells;
    }
  }
  allocateGrid(m_BGrid, num_new_cells, nodes.size());
  
  QVector<vtkIdType> _nodes(grid->GetNumberOfPoints());
  vtkIdType id_new_node = 0;
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    grid->GetPoints()->GetPoint(id_node, x.data());
    m_BGrid->GetPoints()->SetPoint(id_new_node, x.data());
    copyNodeData(grid,id_node,m_BGrid,id_new_node);
    _nodes[id_node] = id_new_node;
    ++id_new_node;
  }
  foreach (vtkIdType id_cell, cells) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = grid->GetCellType(id_cell);
    grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<vtkIdType> new_pts(3);
    if (type_cell == VTK_TRIANGLE) {
      new_pts[0] = _nodes[pts[0]];
      new_pts[1] = _nodes[pts[1]];
      new_pts[2] = _nodes[pts[2]];
      vtkIdType id_new_cell = m_BGrid->InsertNextCell(VTK_TRIANGLE, 3, new_pts.data());
      copyCellData(grid, id_cell, m_BGrid, id_new_cell);
    } else if (type_cell == VTK_QUAD) {
      new_pts[0] = _nodes[pts[0]];
      new_pts[1] = _nodes[pts[1]];
      new_pts[2] = _nodes[pts[2]];
      vtkIdType id_new_cell1 = m_BGrid->InsertNextCell(VTK_TRIANGLE, 3, new_pts.data());
      copyCellData(grid, id_cell, m_BGrid, id_new_cell1);
      new_pts[0] = _nodes[pts[2]];
      new_pts[1] = _nodes[pts[3]];
      new_pts[2] = _nodes[pts[0]];
      vtkIdType id_new_cell2 = m_BGrid->InsertNextCell(VTK_TRIANGLE, 3, new_pts.data());
      copyCellData(grid, id_cell, m_BGrid, id_new_cell2);
    } else {
      EG_BUG;
    }
  }
}

#endif // TRISURFACEPROJECTION_H
