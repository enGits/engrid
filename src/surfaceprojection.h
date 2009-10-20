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

  vec3_t calcGradG(vec3_t x);
  double calcG(vec3_t x);

  vec3_t projectWithLevelSet(vec3_t x);

  bool   projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, const Triangle& T);
  vec3_t projectWithGeometry(vec3_t x, vtkIdType id_node);
  vec3_t correctCurvature(int i_tri, vec3_t r);

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
  vtkIdType addBezierSurface(vtkUnstructuredGrid* bezier, int offset, int N, vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  void writeBezierSurface(vec3_t X_200, vec3_t X_020, vec3_t X_002, vec3_t X_011, vec3_t X_101, vec3_t X_110);
  
  void setupInterpolationGrid();
  
  int getControlPoints_orthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
  int getControlPoints_nonorthogonal(Triangle T, vec3_t& X_011, vec3_t& X_101, vec3_t& X_110);
  
};

template <class C>
void SurfaceProjection::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
{
  setBackgroundGrid_setupGrid(grid, cells);
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
  // create m_Triangles
  m_Triangles.resize(m_BGrid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < m_BGrid->GetNumberOfCells(); ++id_cell) {
    vtkIdType Npts, *pts;
    m_BGrid->GetCellPoints(id_cell, Npts, pts);
    if (Npts == 3) {
      m_BGrid->GetPoints()->GetPoint(pts[0], m_Triangles[id_cell].a.data());
      m_BGrid->GetPoints()->GetPoint(pts[1], m_Triangles[id_cell].b.data());
      m_BGrid->GetPoints()->GetPoint(pts[2], m_Triangles[id_cell].c.data());
      m_Triangles[id_cell].id_a = pts[0];
      m_Triangles[id_cell].id_b = pts[1];
      m_Triangles[id_cell].id_c = pts[2];
      m_Triangles[id_cell].g1 = m_Triangles[id_cell].b - m_Triangles[id_cell].a;
      m_Triangles[id_cell].g2 = m_Triangles[id_cell].c - m_Triangles[id_cell].a;
      m_Triangles[id_cell].g3 = m_Triangles[id_cell].g1.cross(m_Triangles[id_cell].g2);
      m_Triangles[id_cell].A  = 0.5*m_Triangles[id_cell].g3.abs();
      m_Triangles[id_cell].g3.normalise();
      m_Triangles[id_cell].G.column(0, m_Triangles[id_cell].g1);
      m_Triangles[id_cell].G.column(1, m_Triangles[id_cell].g2);
      m_Triangles[id_cell].G.column(2, m_Triangles[id_cell].g3);
      m_Triangles[id_cell].GI = m_Triangles[id_cell].G.inverse();
      m_Triangles[id_cell].smallest_length = (m_Triangles[id_cell].b - m_Triangles[id_cell].a).abs();
      m_Triangles[id_cell].smallest_length = min(m_Triangles[id_cell].smallest_length, (m_Triangles[id_cell].c - m_Triangles[id_cell].b).abs());
      m_Triangles[id_cell].smallest_length = min(m_Triangles[id_cell].smallest_length, (m_Triangles[id_cell].a - m_Triangles[id_cell].c).abs());
    } else {
      EG_ERR_RETURN("only triangles allowed at the moment");
    }
  }

  // compute node normals
  m_NodeNormals.resize(m_BGrid->GetNumberOfPoints());
//   QVector <double> NodeAngles(m_NodeNormals.size());
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node] = vec3_t(0,0,0);
//     NodeAngles[id_node]=0;
  }
  foreach (Triangle T, m_Triangles) {
    double angle_a = GeometryTools::angle(m_BGrid,T.id_c,T.id_a,T.id_b);
    double angle_b = GeometryTools::angle(m_BGrid,T.id_a,T.id_b,T.id_c);
    double angle_c = GeometryTools::angle(m_BGrid,T.id_b,T.id_c,T.id_a);
    double total_angle = angle_a + angle_b + angle_c;
    m_NodeNormals[T.id_a] += angle_a*T.g3;
    m_NodeNormals[T.id_b] += angle_b*T.g3;
    m_NodeNormals[T.id_c] += angle_c*T.g3;
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    m_NodeNormals[id_node].normalise();
  }

  // compute maximum angle per node
  QVector<double> min_cos(m_BGrid->GetNumberOfPoints(), 1.0);
  foreach (Triangle T, m_Triangles) {
    double cosa = T.g3*m_NodeNormals[T.id_a];
    double cosb = T.g3*m_NodeNormals[T.id_b];
    double cosc = T.g3*m_NodeNormals[T.id_c];
    min_cos[T.id_a] = min(cosa, min_cos[T.id_a]);
    min_cos[T.id_b] = min(cosb, min_cos[T.id_b]);
    min_cos[T.id_c] = min(cosc, min_cos[T.id_c]);
  }
  for (vtkIdType id_node = 0; id_node < m_BGrid->GetNumberOfPoints(); ++id_node) {
    double s = sqrt(1.0 - sqr(min(1 - 1e-20, min_cos[id_node])));
    m_EdgeLength[id_node] *= m_RadiusFactor*min_cos[id_node]/s;
  }
}

#endif // SURFACEPROJECTION_H
