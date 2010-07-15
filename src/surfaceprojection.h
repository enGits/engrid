//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "guimainwindow.h"
#include "geometrytools.h"
#include "vtkCharArray.h"
#include "surfaceoperation.h"
#include "surfacealgorithm.h"
#include "triangle.h"
#include "facefinder.h"

class SurfaceProjection : public SurfaceAlgorithm
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
  vtkUnstructuredGrid*      m_FGrid; ///< the foreground grid to project
  QVector<double>           m_EdgeLength;
  QVector<vtkIdType>        m_Cells;
  QVector<vtkIdType>        m_Nodes;
  QVector<vec3_t>           m_NodeNormals; ///< The surface normal at each node of m_BGrid
  QVector<Triangle>         m_Triangles; ///< All triangles of m_BGrid. One for each triangle cell of m_BGrid.
  QVector<double>           m_Radius; ///< Surface radius for mesh resolution.
  QVector<QVector<int> >    m_N2N;
  bool                      m_correctCurvature; ///< Should correctCurvature() be used?
  int                       m_BC;
  bool                      m_RestrictToTriangle;
  double                    m_CritDistance;
  QMap<vtkIdType,vtkIdType> m_Pindex;
  FaceFinder                m_FaceFinder;

protected: // static attributes

  static vtkIdType m_LastPindex;
  virtual vec3_t project(vec3_t x, vtkIdType id_node = -1);

protected: // methods

  virtual void updateBackgroundGridInfo();///< Set up the background grid (triangles, bezier triangles, etc)
  virtual vec3_t correctCurvature(int, vec3_t g_M);
  void searchNewTriangle(vec3_t xp, vtkIdType &id_tri, vec3_t &x_proj, vec3_t &r_proj, bool &on_triangle);
  vtkIdType getProjTriangle(vtkIdType id_node);
  void setProjTriangle(vtkIdType id_node, vtkIdType proj_triangle);
  void computeSurfaceCurvature();

public: // methods

  static long int Nfull;
  static long int Nhalf;

  SurfaceProjection(int bc = 0);
  ~SurfaceProjection();
  
  template <class C> void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells); ///< Set the background grid to use + set it up

  void    setForegroundGrid(vtkUnstructuredGrid* grid);
  virtual vec3_t projectRestricted(vec3_t x, vtkIdType id_node = -1);
  virtual vec3_t projectFree(vec3_t x, vtkIdType id_node = -1);

  vtkUnstructuredGrid* getBGrid() { return m_BGrid; }
  double getRadius(vtkIdType id_node);

  void setCorrectCurvature(bool b) { m_correctCurvature = b; }
  bool getCorrectCurvature()       { return m_correctCurvature; }

public: // static methods

  static void resetPindex() { m_LastPindex = 0; }

};


template <class C>
void SurfaceProjection::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  setBackgroundGrid_setupGrid(grid, cells);
  updateBackgroundGridInfo();
  setForegroundGrid(grid);
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    vtkIdType id_cell = cells[i_cells];
  }
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    m_BGrid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; ++i) {
      setProjTriangle(pts[i], id_cell);
    }
  }
  m_FaceFinder.setMaxNumFaces(100);
  m_FaceFinder.setGrid(m_BGrid);
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
    copyNodeData(grid,id_node,m_BGrid,id_new_node);
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
    vtkIdType id_new_cell = m_BGrid->InsertNextCell(type_cell, N_pts, new_pts);
    copyCellData(grid,id_cell,m_BGrid,id_new_cell);
  }
}

#endif // SURFACEPROJECTION_H
