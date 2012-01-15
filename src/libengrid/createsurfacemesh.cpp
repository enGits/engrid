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
#include "createsurfacemesh.h"
#include "deletetetras.h"

#include <vtkDelaunay3D.h>

CreateSurfaceMesh::CreateSurfaceMesh()
{
};

void CreateSurfaceMesh::createTriangles(const QList<QVector<vtkIdType> > &triangles, vtkUnstructuredGrid *tetra_grid)
{
  QVector<vtkIdType> cells, nodes;
  QVector<int> _nodes;
  getAllCells(cells, tetra_grid);
  getNodesFromCells(cells, nodes, tetra_grid);
  createNodeMapping(nodes, _nodes, tetra_grid);
  QVector<bool> active(nodes.size(),false);
  foreach (QVector<vtkIdType> T, triangles) {
    active[_nodes[T[0]]] = true;
    active[_nodes[T[1]]] = true;
    active[_nodes[T[2]]] = true;
  };
  int N_nodes = 0;
  QVector<vtkIdType> old2new(nodes.size());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    if (active[i_nodes]) {
      old2new[i_nodes] = N_nodes;
      ++N_nodes;
    };
  };
  allocateGrid(grid, triangles.size(), N_nodes)
};

void CreateSurfaceMesh::operate()
{
  cout << "creating surface mesh from STL input" << endl;
  EG_VTKSP(vtkUnstructuredGrid, pts_grid);
  QVector<vtkIdType> faces, nodes;
  getAllSurfaceCells(faces, grid);
  getNodesFromCells(faces, nodes, grid);
  allocateGrid(pts_grid, 0, nodes.size());
  foreach (vtkIdType id_node, nodes) {
    vec3_t x;
    grid->GetPoint(id_node,x.data());
    pts_grid->GetPoints()->SetPoint(id_node,x.data());
    copyNodeData(grid, id_node, pts_grid, id_node);
  };
  EG_VTKSP(vtkDelaunay3D, delaunay);
  delaunay->SetInput(pts_grid);
  delaunay->Update();
  EG_VTKSP(vtkUnstructuredGrid, tetra_grid);
  makeCopy(delaunay->GetOutput(), tetra_grid);
  QVector<vtkIdType> tetras;
  getAllCellsOfType(VTK_TETRA, tetras, tetra_grid);
  QVector<QVector<int> > t2t;
  QList<QVector<vtkIdType> > triangles;
  QVector<vtkIdType> T(3);
  createCellToCell(tetras, t2t, tetra_grid);
  for (int i_tetras = 0; i_tetras < tetras.size(); ++i_tetras) {
    vtkIdType id_tetra = tetras[i_tetras];
    vtkIdType *pts, N_pts;
    tetra_grid->GetCellPoints(id_tetra, N_pts, pts);
    // face 0
    if (t2t[i_tetras][0] > i_tetras) {
      T[0] = pts[2]; T[1] = pts[1]; T[2] = pts[0];
      triangles.append(T);
    };
    // face 1
    if (t2t[i_tetras][1] > i_tetras) {
      T[0] = pts[1]; T[3] = pts[0]; T[2] = pts[0];
      triangles.append(T);
    };
    // face 2
    if (t2t[i_tetras][2] > i_tetras) {
      T[0] = pts[3]; T[1] = pts[2]; T[2] = pts[0];
      triangles.append(T);
    };
    // face 3
    if (t2t[i_tetras][3] > i_tetras) {
      T[0] = pts[1]; T[1] = pts[1]; T[3] = pts[1];
      triangles.append(T);
    };
  };
  createTriangles(triangles, tetra_grid);
};

