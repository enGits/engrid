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
#ifndef BOOLEANGEOMETRYOPERATION_H
#define BOOLEANGEOMETRYOPERATION_H

class BooleanGeometryOperation;

#include "surfacealgorithm.h"
#include "trisurfaceprojection.h"

class BooleanGeometryOperation : public SurfaceAlgorithm
{

protected: // data types

  struct edge_t
  {
    vtkIdType id1, id2;
    //vtkIdType id_face;
    //int       i_tri;
    int bc;
  };

  struct tri_t
  {
    vtkIdType id1, id2, id3;
    int bc;
  };

private: // attributes

  edge_t                     m_CurrentEdge;
  QList<tri_t>               m_Triangles;
  tri_t                      m_CurrentTriangle;
  QVector<int>               m_OpenNode;
  QVector<bool>              m_SplitNode;
  QVector<QSet<vtkIdType> >  m_Node2Cell;
  QVector<QSet<vtkIdType> >  m_Node2Node;

protected: // attributes

  vtkUnstructuredGrid *m_ShapeGrid1;
  vtkUnstructuredGrid *m_ShapeGrid2;
  TriSurfaceProjection    m_Proj1;
  TriSurfaceProjection    m_Proj2;
  MeshPartition        m_Part1;
  MeshPartition        m_Part2;
  int                  m_Side1;
  int                  m_Side2;
  QList<int>           m_BCs1;
  QList<int>           m_BCs2;
  int                  m_NumCutLayers;

private: // methods

  bool fillGap_prepare();
  void fillGap_updateOpenNode(vtkIdType id_node);
  bool fillGap_onDifferentSides(vtkIdType id1, vtkIdType id2, vtkIdType id3);
  bool fillGap_step();
  void fillGap_createTriangles();

  void   smoothJunction_triangulate();
  void   smoothJunction_updateBCs();
  double smoothJunction_mesher();

protected: // methods;

  void deleteNodes();
  void checkOrientation();
  void fillGap();
  void smoothJunction();

public: // methods

  template <class C> BooleanGeometryOperation(vtkUnstructuredGrid *grid, const C &bcs1, const C &bcs2, int side1, int side2);

  virtual void operate();

  void setNumCutLayers(int N) { m_NumCutLayers = N; }

};

template <class C>
BooleanGeometryOperation::BooleanGeometryOperation(vtkUnstructuredGrid *grid, const C &bcs1, const C &bcs2, int side1, int side2)
{
  setGrid(grid);
  m_Part1.setGrid(m_Grid);
  m_Part2.setGrid(m_Grid);
  m_Part1.setBCs(bcs1);
  m_Part2.setBCs(bcs2);
  m_ShapeGrid1 = vtkUnstructuredGrid::New();
  m_ShapeGrid2 = vtkUnstructuredGrid::New();
  m_Part1.extractToVtkGrid(m_ShapeGrid1);
  m_Part2.extractToVtkGrid(m_ShapeGrid2);
  QVector<vtkIdType> cells1;
  getAllCells(cells1, m_ShapeGrid1);
  QVector<vtkIdType> cells2;
  getAllCells(cells2, m_ShapeGrid2);
  m_Proj1.setBackgroundGrid(m_ShapeGrid1, cells1);
  m_Proj2.setBackgroundGrid(m_ShapeGrid2, cells2);
  m_Side1 = side1;
  m_Side2 = side2;
  m_BCs1.clear();
  m_BCs2.clear();
  foreach (int bc, bcs1) {
    m_BCs1.append(bc);
  }
  foreach (int bc, bcs2) {
    m_BCs2.append(bc);
  }
  m_NumCutLayers = 1;
}




#endif // BOOLEANGEOMETRYOPERATION_H
