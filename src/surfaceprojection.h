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
#include "geometrytools.h"
#include "vtkCharArray.h"
#include "surfaceoperation.h"
#include "surfacealgorithm.h"

class SurfaceProjection : public SurfaceAlgorithm
{

private: // data-types

  struct Triangle
  {
    vtkIdType id_a, id_b, id_c;
    vec3_t a, b, c;
    vec3_t g1, g2, g3;
    mat3_t G, GI;
    double A;
    double smallest_length;
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
  
  vtkUnstructuredGrid*   m_InterpolationGrid;
  vtkUnstructuredGrid*   m_BezierGrid;
  
  QVector<vtkIdType>     m_ProjTriangles;
  vtkUnstructuredGrid*   m_FGrid;
  QVector<double>        m_EdgeLength;
  QVector<vtkIdType>     m_Cells;
  QVector<vtkIdType>     m_Nodes;
  QVector<vec3_t>        m_NodeNormals;
  QVector<Triangle>      m_Triangles;
  QVector<QVector<int> > m_N2N;
  Octree                 m_OTGrid;
  QVector<double>        m_G;
  QVector<bool>          m_GSet;
  double                 m_Relax;
  double                 m_DistExp;
  double                 m_DistWeight;
  double                 m_DirWeight;
  double                 m_DirExp;
  double                 m_WeightOffset;
  double                 m_MinOTLength;
  int                    m_MaxOTCells;
  double                 m_Length;
  int                    m_MaxIter;
  double                 m_ConvLimit;
  double                 m_RadiusFactor;
  bool                   m_UseLevelSet;
  int                    m_NumDirect;
  int                    m_NumFull;

private: // methods

  template <class C>
  void setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells);
  void setBackgroundGrid_initOctree();
  void setBackgroundGrid_refineFromNodes();
  void setBackgroundGrid_refineFromEdges();
  void setBackgroundGrid_refineFromFaces();
  void setBackgroundGrid_initLevelSet();
  void setBackgroundGrid_computeLevelSet();

  void updateBackgroundGridInfo();
  
  vec3_t calcGradG(vec3_t x);
  double calcG(vec3_t x);

  vec3_t projectWithLevelSet(vec3_t x);

  bool   projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, const Triangle& T);
  vec3_t projectWithGeometry(vec3_t x, vtkIdType id_node);
  vec3_t correctCurvature(int i_tri, vec3_t r);
  vec3_t cylinder(vec3_t center, double radius, vec3_t g_M);
  vec3_t cylinder(vec3_t center, double radius, int i_tri, vec3_t r);
  
public: // methods

  SurfaceProjection();

  template <class C>
  void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells);

  void setForegroundGrid(vtkUnstructuredGrid* grid) {m_FGrid = grid; }

  vec3_t project(vec3_t x, vtkIdType id_node = -1);
  int getNumOctreeCells() { return m_OTGrid.getNumCells(); }
  void writeOctree(QString file_name);
  bool usesLevelSet() { return m_UseLevelSet; };

  int getNumDirectProjections() { return m_NumDirect; }
  int getNumFullSearches() { return m_NumFull; }

  void writeGridWithNormals();
//   vtkIdType addBezierSurface(vtkUnstructuredGrid* bezier, int offset, int N, vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  void writeBezierSurface(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  
  void setupInterpolationGrid();
  
  int getControlPoints_orthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
  int getControlPoints_nonorthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
  
//   mat3_t Jacobian_Matrix(double x, double y);
  
};

template <class C>
void SurfaceProjection::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  setBackgroundGrid_setupGrid(grid, cells);
  updateBackgroundGridInfo();
  if (m_UseLevelSet) {
    setBackgroundGrid_initOctree();
    setBackgroundGrid_refineFromNodes();
    setBackgroundGrid_refineFromEdges();
    setBackgroundGrid_refineFromFaces();
    setBackgroundGrid_computeLevelSet();
  }
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
    m_BGrid->InsertNextCell(type_cell, N_pts, new_pts);
  }
}

#endif // SURFACEPROJECTION_H
