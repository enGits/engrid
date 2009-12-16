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
#include "triangle.h"
#include "beziertriangle.h"

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

private: // attributes

  vtkUnstructuredGrid*   m_BGrid; ///< background grid used for projection and interpolation of the bezier surface
  vtkUnstructuredGrid*   m_InterpolationGrid;
  vtkUnstructuredGrid*   m_BezierGrid;
  
public:
  /// get m_BGrid
  vtkUnstructuredGrid* getBGrid() { return m_BGrid; }
  /// get m_InterpolationGrid
  vtkUnstructuredGrid* getInterpolationGrid() { return m_InterpolationGrid; }
    /// get m_BezierGrid
  vtkUnstructuredGrid* getBezierGrid() { return m_BezierGrid; }
  
private:
  /// A vector associating each node of m_FGrid with a Triangle (index for m_Triangles) on which it should be projected (closest triangle from m_Triangles).
  QVector<vtkIdType>     m_ProjTriangles;
  
  vtkUnstructuredGrid*    m_FGrid; ///< The foreground grid to project.
  QVector<double>         m_EdgeLength;
  QVector<vtkIdType>      m_Cells;
  QVector<vtkIdType>      m_Nodes;
  QVector<vec3_t>         m_NodeNormals; ///< The surface normal at each node of m_BGrid
  QMap < pair <vtkIdType, vtkIdType>, vec3_t > m_ControlPoints;
  QVector<Triangle>       m_Triangles; ///< All triangles of m_BGrid. One for each triangle cell of m_BGrid.
  QVector<BezierTriangle> m_BezierTriangles; ///< The bezier triangle corresponding to m_Triangles
  QVector<QVector<int> >  m_N2N;
  Octree                  m_OTGrid;
  QVector<double>         m_G;
  QVector<bool>           m_GSet;
  double                  m_Relax;
  double                  m_DistExp;
  double                  m_DistWeight;
  double                  m_DirWeight;
  double                  m_DirExp;
  double                  m_WeightOffset;
  double                  m_MinOTLength;
  int                     m_MaxOTCells;
  double                  m_Length;
  int                     m_MaxIter;
  double                  m_ConvLimit;
  double                  m_RadiusFactor;
  bool                    m_UseLevelSet;
  int                     m_NumDirect;
  int                     m_NumFull;

  bool m_correctCurvature; ///< Should correctCurvature() be used?

public:
  void setCorrectCurvature(bool b) { m_correctCurvature = b; }
  bool getCorrectCurvature() { return m_correctCurvature; }
  
// variables for exact projection surfaces
public:
  int m_ExactMode;
  vec3_t m_center;
  vec3_t m_Rx;
  vec3_t m_Ry;
  vec3_t m_Rz;
  vec3_t cylinder(vec3_t center, double radius, vec3_t g_M);
  vec3_t cylinder(vec3_t center, double radius, int i_tri, vec3_t r);
  vec3_t ellipsoid(vec3_t M);
  vec3_t ellipse(vec3_t M);
  vec3_t rectangle(vec3_t M);
  vec3_t cuboid(vec3_t M);
  vec3_t cylinder(vec3_t M);
  
private: // methods

  template <class C>
  void setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells); ///< copy the cells from grid to m_BGrid
  void setBackgroundGrid_initOctree();
  void setBackgroundGrid_refineFromNodes();
  void setBackgroundGrid_refineFromEdges();
  void setBackgroundGrid_refineFromFaces();
  void setBackgroundGrid_initLevelSet();
  void setBackgroundGrid_computeLevelSet();

  void updateBackgroundGridInfo();///< Set up the background grid (triangles, bezier triangles, etc)
  void updateBackgroundGridInfo_original();
  
  vec3_t calcGradG(vec3_t x);
  double calcG(vec3_t x);

  vec3_t projectWithLevelSet(vec3_t x);

  vec3_t projectWithGeometry(vec3_t x, vtkIdType id_node); ///< project onto m_BGrid, eventually using correctCurvature
  vec3_t projectWithGeometry_original(vec3_t x, vtkIdType id_node);
    
  vec3_t correctCurvature1(int i_tri, vec3_t g_M); ///< correct curvature by using double interpolation
  vec3_t correctCurvature2(int i_tri, vec3_t g_M); ///< correct curvature by using bezier surfaces
  
  vec3_t getEdgeNormal(vtkIdType id_node1, vtkIdType id_node2);
  
public: // methods

  SurfaceProjection(); ///< Constructor
  ~SurfaceProjection(); ///< Destructor
  
  template <class C>
  void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells); ///< Set the background grid to use + set it up

  void setForegroundGrid(vtkUnstructuredGrid* grid) {
    m_FGrid = grid;
    m_ProjTriangles.clear(); // this makes sure a full search is run everytime a new foreground grid is set.
  } ///< set m_FGrid

  vec3_t project(vec3_t x, vtkIdType id_node = -1);
  int getNumOctreeCells() { return m_OTGrid.getNumCells(); }
  void writeOctree(QString file_name);
  bool usesLevelSet() { return m_UseLevelSet; }; ///< Set whether or not to use the level set method

  int getNumDirectProjections() { return m_NumDirect; }
  int getNumFullSearches() { return m_NumFull; }

  void writeGridWithNormals(QString filename);
  void writeInterpolationGrid(QString filename);
  void writeTriangleGrid(QString filename);
  
  int getControlPoints_orthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110, double Lmax); ///< get the orthogonal control points
  int getControlPoints_nonorthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110, double Lmax); ///< get the non-orthogonal control points
  int limitControlPoints(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
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
    vtkIdType id_new_cell = m_BGrid->InsertNextCell(type_cell, N_pts, new_pts);
    copyCellData(grid,id_cell,m_BGrid,id_new_cell);
  }
}

#endif // SURFACEPROJECTION_H
