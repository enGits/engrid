//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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

void CreateSurfaceMesh::operate()
{
  cout << "creating surface mesh from STL input" << endl;
  EG_VTKSP(vtkUnstructuredGrid, pts_grid);
  QVector<vtkIdType> faces, nodes;
  getAllSurfaceCells(faces, grid);
  getNodesFromCells(faces, nodes, grid);
  allocateGrid(pts_grid, faces.size(), nodes.size());
  EG_VTKSP(vtkDelaunay3D, delaunay);
  delaunay->SetInput(pts_grid);
  delaunay->Update();
  makeCopy(delaunay->GetOutput(), grid);
};

