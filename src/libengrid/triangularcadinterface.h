//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                     +
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
#ifndef TRIANGULARCADINTERFACE_H
#define TRIANGULARCADINTERFACE_H

#include "cadinterface.h"
#include "triangle.h"
#include "facefinder.h"
#include "surfacealgorithm.h"

class TriangularCadInterface : public CadInterface, public SurfaceAlgorithm
{

private: // methods

  template <class C>
  void setBackgroundGrid_setupGrid(vtkUnstructuredGrid* grid, const C& cells); ///< copy the cells from grid to m_BGrid

  void updateBackgroundGridInfo();      ///< Set up the background grid (triangles, bezier triangles, etc)
  void computeSurfaceCurvature();

protected: // attributes

  vtkUnstructuredGrid*      m_BGrid;       ///< the background grid defining the geometry
  MeshPartition             m_BPart;
  QVector<vtkIdType>        m_Cells;
  QVector<vtkIdType>        m_Nodes;
  QVector<vec3_t>           m_NodeNormals; ///< The surface normal at each node of m_BGrid
  QVector<Triangle>         m_Triangles;   ///< All triangles of m_BGrid. One for each triangle cell of m_BGrid.
  QVector<double>           m_Radius;      ///< Surface radius for mesh resolution.
  QVector<QVector<int> >    m_N2N;
  double                    m_CritDistance;
  FaceFinder                m_FaceFinder;


public:

  TriangularCadInterface();

  virtual HitType      shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r);
  virtual PositionType position(vec3_t x, vec3_t n);

  template <class C> void setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells); ///< Set the background grid to use + set it up

};

template <class C>
void TriangularCadInterface::setBackgroundGrid(vtkUnstructuredGrid* grid, const C& cells)
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
  updateBackgroundGridInfo();
  m_FaceFinder.setMaxNumFaces(10);
  m_FaceFinder.setGrid(m_BGrid);
}

#endif // TRIANGULARCADINTERFACE_H
