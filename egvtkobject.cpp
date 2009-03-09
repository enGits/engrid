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
#include "egvtkobject.h"
#include "guimainwindow.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellLinks.h>
#include <vtkCellType.h>
#include <vtkIdList.h>
#include <vtkCell.h>

void EgVtkObject::computeNormals
(
  QVector<vec3_t>     &cell_normals,
  QVector<vec3_t>     &node_normals,
  QVector<vtkIdType>  &cells,
  QVector<vtkIdType>  &nodes,
  vtkUnstructuredGrid *grid
)
{
  using namespace GeometryTools;
  
  cell_normals.resize(cells.size());
  node_normals.fill(vec3_t(0,0,0), nodes.size());
  QVector<int> g2s;
  createNodeMapping(nodes, g2s, grid);
  for (int i_cell = 0; i_cell < cells.count(); ++i_cell) {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType *pts;
    vtkIdType npts;
    grid->GetCellPoints(id_cell, npts, pts);
    cell_normals[i_cell] = cellNormal(grid, id_cell);
    cell_normals[i_cell].normalise();
    for (int i_pts = 0; i_pts < npts; ++i_pts) {
      if (g2s[pts[i_pts]] != -1) {
        node_normals[g2s[pts[i_pts]]] += cell_normals[i_cell];
      };
    };
  };
  for (int i_node = 0; i_node < nodes.count(); ++i_node) {
    node_normals[i_node].normalise();
    //cout << node_normals[i_node] << endl;
  };
};

void EgVtkObject::createNodeMapping
(
  QVector<vtkIdType>  &nodes,
  QVector<int>        &_nodes,
  vtkUnstructuredGrid *grid
)
{
  _nodes.fill(-1,grid->GetNumberOfPoints());
  vtkIdType id;
  int i_sub = 0;
  foreach(id,nodes) {
    _nodes[id] = i_sub;
    ++i_sub;
  };
};

void EgVtkObject::createCellMapping
(
  QVector<vtkIdType>  &cells,
  QVector<int>        &_cells,
  vtkUnstructuredGrid *grid
)
{
  _cells.fill(-1,grid->GetNumberOfCells());
  vtkIdType id;
  int i_sub = 0;
  foreach(id,cells) {
    _cells[id] = i_sub;
    ++i_sub;
  };
};

void EgVtkObject::createNodeToBcMapping
(
  QVector<QSet<int> > &bcs,
  vtkUnstructuredGrid *grid
)
{
  bcs.fill(QSet<int>(), grid->GetNumberOfPoints());
  grid->BuildLinks();
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType nodeId = 0; nodeId < grid->GetNumberOfPoints(); ++nodeId) {
    int Ncells = grid->GetCellLinks()->GetNcells(nodeId);
    for (int i = 0; i < Ncells; ++i) {
      vtkIdType id_cell = grid->GetCellLinks()->GetCells(nodeId)[i];
      vtkIdType ct = grid->GetCellType(id_cell);
      if ((ct == VTK_TRIANGLE) || (ct = VTK_QUAD)) {
        if (cell_code->GetValue(id_cell) > 0) {
          bcs[nodeId].insert(cell_code->GetValue(id_cell));
        };
      };
    };
  };
};

void EgVtkObject::createNodeToCell
(
  QVector<vtkIdType>  &cells,
  QVector<vtkIdType>  &nodes,
  QVector<int>        &_nodes,
  QVector<QSet<int> > &n2c,
  vtkUnstructuredGrid *grid
)
{
  n2c.fill(QSet<int>(), nodes.size());
  for (vtkIdType i_cells = 0; i_cells < cells.size(); ++i_cells) {
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(cells[i_cells], Npts, pts);
    for (int i_pts = 0; i_pts < Npts; ++i_pts) {
      n2c[_nodes[pts[i_pts]]].insert(i_cells);
    };
  };
};

void EgVtkObject::addToN2N
(
  QVector<QSet<int> > &n2n, 
  int                  n1, 
  int                  n2
)
{
  n2n[n1].insert(n2);
  n2n[n2].insert(n1);
};

void EgVtkObject::createNodeToNode
(
  QVector<vtkIdType>  &cells,
  QVector<vtkIdType>  &nodes,
  QVector<int>        &_nodes,
  QVector<QSet<int> > &n2n,
  vtkUnstructuredGrid *grid
)
{
  n2n.fill(QSet<int>(), nodes.size());
  foreach (vtkIdType id_cell, cells) {
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(id_cell, Npts, pts);
    vector<int> n(Npts);
    for (int i = 0; i < Npts; ++i) {
      n[i] = _nodes[pts[i]];
    };
    vtkIdType cellType = grid->GetCellType(id_cell);
    if (cellType == VTK_TRIANGLE) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[2], n[0]);
    } else if (cellType == VTK_QUAD) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[2], n[3]);
      addToN2N(n2n, n[3], n[0]);
    } else if (cellType == VTK_TETRA) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[2]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[3]);
      addToN2N(n2n, n[2], n[3]);
    } else if (cellType == VTK_PYRAMID) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[0], n[4]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[4]);
      addToN2N(n2n, n[2], n[3]);
      addToN2N(n2n, n[2], n[4]);
      addToN2N(n2n, n[3], n[4]);
    } else if (cellType == VTK_WEDGE) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[2]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[4]);
      addToN2N(n2n, n[2], n[5]);
      addToN2N(n2n, n[3], n[4]);
      addToN2N(n2n, n[3], n[5]);
      addToN2N(n2n, n[4], n[5]);
    } else if (cellType == VTK_HEXAHEDRON) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[0], n[4]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[5]);
      addToN2N(n2n, n[2], n[3]);
      addToN2N(n2n, n[2], n[6]);
      addToN2N(n2n, n[3], n[7]);
      addToN2N(n2n, n[4], n[5]);
      addToN2N(n2n, n[4], n[7]);
      addToN2N(n2n, n[5], n[6]);
      addToN2N(n2n, n[6], n[7]);
    };
  };
};

void EgVtkObject::getAllCells
(
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    ++N;
  };
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cells[N] = id_cell;
    ++N;
  };
};

void EgVtkObject::getAllCellsOfType
(
  vtkIdType            type,
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (grid->GetCellType(id_cell) == type) {
      ++N;
    };
  };
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (grid->GetCellType(id_cell) == type) {
      cells[N] = id_cell;
      ++N;
    };
  };
};


void EgVtkObject::getAllVolumeCells
(
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      ++N;
    };
  };
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      cells[N] = id_cell;
      ++N;
    };
  };
};

void EgVtkObject::getAllSurfaceCells
(
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      ++N;
    };
  };
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      cells[N] = id_cell;
      ++N;
    };
  };
};

void EgVtkObject::getSurfaceCells
(
  QSet<int>           &bcs,
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      if (bcs.contains(cell_code->GetValue(id_cell))) {
        ++N;
      };
    };
  };
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      if (bcs.contains(cell_code->GetValue(id_cell))) {
        cells[N] = id_cell;
        ++N;
      };
    };
  };
};

void EgVtkObject::getSurfaceNodes
(
  QSet<int>           &bcs,
  QSet <vtkIdType> &SelectedNodes,
  vtkUnstructuredGrid *grid
)
{
  QVector<vtkIdType> SelectedCells;
  getSurfaceCells(bcs, SelectedCells, grid);
  
  SelectedNodes.clear();
  foreach(vtkIdType id_cell, SelectedCells)
  {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for(int i=0;i<N_pts;i++)
    {
      SelectedNodes.insert(pts[i]);
    }
  }
};

void EgVtkObject::addToC2C
(
  vtkIdType               id_cell,
  QVector<int>           &_cells,
  QVector<QVector<int> > &c2c,
  int                     j,
  vtkIdList              *nds,
  vtkIdList              *cls,
  vtkUnstructuredGrid    *grid
)
{
  c2c[_cells[id_cell]][j] = -1;
  grid->GetCellNeighbors(id_cell, nds, cls);
  for (int i = 0; i < cls->GetNumberOfIds(); ++i) {
    if (cls->GetId(i) != id_cell) {
      if (_cells[cls->GetId(i)] != -1) {
        c2c[_cells[id_cell]][j] = _cells[cls->GetId(i)];
      };
    };
  };
};


void EgVtkObject::createCellToCell
(
  QVector<vtkIdType>     &cells,
  QVector<QVector<int> > &c2c,
  vtkUnstructuredGrid    *grid
)
{
  // GetCellNeighbors(vtkIdType id_cell, vtkIdList *ptIds, vtkIdList *id_cells)
  grid->BuildLinks();
  QVector<int> _cells;
  createCellMapping(cells, _cells, grid);
  c2c.fill(QVector<int>(), cells.size());
  EG_VTKSP(vtkIdList, nds);
  EG_VTKSP(vtkIdList, cls);
  for (int i = 0; i < cells.size(); ++i) {
    vtkIdType id_cell = cells[i];
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(id_cell, Npts, pts);
    if (grid->GetCellType(id_cell) == VTK_TRIANGLE) {
      c2c[i].resize(3);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[0]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
    } else if (grid->GetCellType(id_cell) == VTK_QUAD) {
      c2c[i].resize(4);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[0]);
      addToC2C(id_cell, _cells, c2c, 3, nds, cls, grid);
    } else if (grid->GetCellType(id_cell) == VTK_TETRA) {
      c2c[i].resize(4);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 3, nds, cls, grid);
    } else if (grid->GetCellType(id_cell) == VTK_PYRAMID) {
      c2c[i].resize(5);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[4]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[4]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[4]);
      addToC2C(id_cell, _cells, c2c, 3, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[4]);
      addToC2C(id_cell, _cells, c2c, 4, nds, cls, grid);
    } else if (grid->GetCellType(id_cell) == VTK_WEDGE) {
      c2c[i].resize(5);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[4]);
      nds->InsertNextId(pts[5]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[4]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[4]);
      nds->InsertNextId(pts[5]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 3, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[5]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 4, nds, cls, grid);
    } else if (grid->GetCellType(id_cell) == VTK_HEXAHEDRON) {
      c2c[i].resize(6);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[1]);
      addToC2C(id_cell, _cells, c2c, 0, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[4]);
      nds->InsertNextId(pts[5]);
      nds->InsertNextId(pts[6]);
      nds->InsertNextId(pts[7]);
      addToC2C(id_cell, _cells, c2c, 1, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[5]);
      nds->InsertNextId(pts[4]);
      addToC2C(id_cell, _cells, c2c, 2, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[3]);
      nds->InsertNextId(pts[7]);
      nds->InsertNextId(pts[6]);
      nds->InsertNextId(pts[2]);
      addToC2C(id_cell, _cells, c2c, 3, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[0]);
      nds->InsertNextId(pts[4]);
      nds->InsertNextId(pts[7]);
      nds->InsertNextId(pts[3]);
      addToC2C(id_cell, _cells, c2c, 4, nds, cls, grid);
      nds->Reset();
      nds->InsertNextId(pts[1]);
      nds->InsertNextId(pts[2]);
      nds->InsertNextId(pts[6]);
      nds->InsertNextId(pts[5]);
      addToC2C(id_cell, _cells, c2c, 5, nds, cls, grid);
    };
  };
};

void EgVtkObject::getNodesFromCells
(
  QVector<vtkIdType>  &cells,
  QVector<vtkIdType>  &nodes,
  vtkUnstructuredGrid *grid
)
{
  QSet<vtkIdType> ex_nodes;
  vtkIdType id_cell;
  foreach(id_cell, cells) {
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(id_cell, Npts, pts);
    for (int i = 0; i < Npts; ++i) {
      ex_nodes.insert(pts[i]);
    };
  };
  nodes.resize(ex_nodes.size());
  {
    int j = 0;
    vtkIdType i;
    foreach(i,ex_nodes) {
      nodes[j] = i;
      ++j;
    };
  };
};

bool EgVtkObject::isVolume(vtkUnstructuredGrid *grid, vtkIdType id_cell)
{
  bool isVol = false;
  if      (grid->GetCellType(id_cell) == VTK_TETRA)      isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_PYRAMID)    isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_WEDGE)      isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_HEXAHEDRON) isVol = true;
  return isVol;
};

bool EgVtkObject::isSurface(vtkUnstructuredGrid *grid, vtkIdType id_cell)
{
  bool isSurf = false;
  if      (grid->GetCellType(id_cell) == VTK_TRIANGLE) isSurf = true;
  else if (grid->GetCellType(id_cell) == VTK_QUAD)     isSurf = true;
  return isSurf;
};

void EgVtkObject::UpdateCellIndex(vtkUnstructuredGrid *grid)
{
  if (!grid->GetCellData()->GetArray("cell_index")) {
    EG_VTKSP(vtkLongArray_t, cell_index);
    cell_index->SetName("cell_index");
    cell_index->SetNumberOfValues(grid->GetNumberOfCells());
    grid->GetCellData()->AddArray(cell_index);
  };
  EG_VTKDCC(vtkLongArray_t, cell_index, grid, "cell_index");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cell_index->SetValue(id_cell, id_cell);
  };
};

void EgVtkObject::UpdateNodeIndex(vtkUnstructuredGrid *grid)
{
  if (!grid->GetPointData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, node_index);
    node_index->SetName("node_index");
    node_index->SetNumberOfValues(grid->GetNumberOfPoints());
    grid->GetPointData()->AddArray(node_index);
  };
  EG_VTKDCN(vtkLongArray_t, node_index, grid, "node_index");
  for (vtkIdType pointId = 0; pointId < grid->GetNumberOfPoints(); ++pointId) {
    node_index->SetValue(pointId, pointId);
  };
};

void EgVtkObject::addToPolyData
(
  QVector<vtkIdType>  &cells,
  vtkPolyData         *pdata,
  vtkUnstructuredGrid *grid
)
{
  UpdateCellIndex(grid);
  UpdateNodeIndex(grid);
  QVector<vtkIdType> nodes;
  QVector<int>       _nodes;
  getNodesFromCells(cells, nodes, grid);
  createNodeMapping(nodes, _nodes, grid);
  EG_VTKSP(vtkDoubleArray, pcoords);
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfTuples(nodes.size());
  EG_VTKSP(vtkPoints, points);
  points->SetData(pcoords);
  pdata->SetPoints(points);
  pdata->Allocate(cells.size());
  if (!pdata->GetCellData()->GetArray("cell_index")) {
    EG_VTKSP(vtkLongArray_t, cell_index);
    cell_index->SetName("cell_index");
    //cell_index->SetNumberOfValues(cells.size());
    pdata->GetCellData()->AddArray(cell_index);
  };
  if (!pdata->GetPointData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, node_index);
    node_index->SetName("node_index");
    //node_index->SetNumberOfValues(nodes.size());
    pdata->GetPointData()->AddArray(node_index);
  };
  EG_VTKDCC(vtkLongArray_t, pd_cell_index, pdata, "cell_index");
  EG_VTKDCN(vtkLongArray_t, pd_node_index, pdata, "node_index");
  pd_cell_index->SetNumberOfValues(cells.size());
  pd_node_index->SetNumberOfValues(nodes.size());
  for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType cellType = grid->GetCellType(id_cell);
    if ((cellType != VTK_TRIANGLE) && (cellType != VTK_QUAD)) {
      EG_ERR_RETURN("unsupported cell type for this operation");
    };
    vtkIdType Npts, *pts;
    grid->GetCellPoints(id_cell, Npts, pts);
    vtkIdType *new_pts = new vtkIdType[Npts];
    for (int i = 0; i < Npts; ++i) {
      new_pts[i] = _nodes[pts[i]];
    };
    vtkIdType newCellId = pdata->InsertNextCell(cellType, Npts, new_pts);
    pd_cell_index->SetValue(newCellId, id_cell);
    delete [] new_pts;
  };
  for (int i_node = 0; i_node < nodes.size(); ++i_node) {
    vec3_t x;
    grid->GetPoints()->GetPoint(nodes[i_node], x.data());
    pdata->GetPoints()->SetPoint(i_node, x.data());
    pd_node_index->SetValue(i_node, nodes[i_node]);
  };
};

#define EGVTKOBJECT_COPYCELLDATA(FIELD,TYPE) \
{ \
if (old_grid->GetCellData()->GetArray(FIELD)) { \
EG_VTKDCC(TYPE, var1, old_grid, FIELD); \
EG_VTKDCC(TYPE, var2, new_grid, FIELD); \
var2->SetValue(newId, var1->GetValue(oldId)); \
}; \
};

void EgVtkObject::copyCellData
(
  vtkUnstructuredGrid *old_grid, 
  vtkIdType            oldId, 
  vtkUnstructuredGrid *new_grid,
  vtkIdType            newId
)
{
  EGVTKOBJECT_COPYCELLDATA("vtk_type",      vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_code",     vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_index",    vtkLongArray_t);
  EGVTKOBJECT_COPYCELLDATA("cell_err_tet",  vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_pria", vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_prib", vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_pric", vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_prid", vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_prie", vtkDoubleArray);
  EGVTKOBJECT_COPYCELLDATA("cell_err_prif", vtkDoubleArray);
};

#define EGVTKOBJECT_COPYNODEDATA(FIELD,TYPE) \
{ \
if (old_grid->GetPointData()->GetArray(FIELD)) { \
EG_VTKDCN(TYPE, var1, old_grid, FIELD); \
EG_VTKDCN(TYPE, var2, new_grid, FIELD); \
var2->SetValue(newId, var1->GetValue(oldId)); \
}; \
};

void EgVtkObject::copyNodeData
(
  vtkUnstructuredGrid *old_grid, 
  vtkIdType            oldId, 
  vtkUnstructuredGrid *new_grid,
  vtkIdType            newId
)
{
  EGVTKOBJECT_COPYNODEDATA("node_status", vtkIntArray);
  EGVTKOBJECT_COPYNODEDATA("node_layer",  vtkIntArray);
  EGVTKOBJECT_COPYNODEDATA("node_index",  vtkLongArray_t);
};

#define EGVTKOBJECT_CREATECELLFIELD(FIELD,TYPE,OW) \
if (!grid->GetCellData()->GetArray(FIELD)) { \
EG_VTKSP(TYPE, var); \
var->SetName(FIELD); \
var->SetNumberOfValues(Ncells); \
grid->GetCellData()->AddArray(var); \
} else if (OW) { \
EG_VTKDCC(TYPE, var, grid, FIELD); \
var->SetNumberOfValues(Ncells); \
};

#define EGVTKOBJECT_CREATENODEFIELD(FIELD,TYPE,OW) \
if (!grid->GetPointData()->GetArray(FIELD)) { \
EG_VTKSP(TYPE, var); \
var->SetName(FIELD); \
var->SetNumberOfValues(Nnodes); \
grid->GetPointData()->AddArray(var); \
for (int i = 0; i < grid->GetNumberOfPoints(); ++i) var->SetValue(i,0); \
} else if (OW) { \
EG_VTKDCN(TYPE, var, grid, FIELD); \
var->SetNumberOfValues(Nnodes); \
for (int i = 0; i < grid->GetNumberOfPoints(); ++i) var->SetValue(i,0); \
};

void EgVtkObject::createBasicFields
(
  vtkUnstructuredGrid *grid,
  vtkIdType            Ncells, 
  vtkIdType            Nnodes,
  bool                 overwrite
)
{
  createBasicNodeFields(grid, Nnodes, overwrite);
  createBasicCellFields(grid, Ncells, overwrite);
};

void EgVtkObject::createBasicCellFields
(
  vtkUnstructuredGrid *grid,
  vtkIdType            Ncells, 
  bool                 overwrite
)
{
  EGVTKOBJECT_CREATECELLFIELD("vtk_type" ,     vtkIntArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_code",     vtkIntArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_index",    vtkLongArray_t, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_tet",  vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_pria", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_prib", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_pric", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_prid", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_prie", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_err_prif", vtkDoubleArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_VA",       vtkDoubleArray, overwrite);
};

void EgVtkObject::createBasicNodeFields
(
  vtkUnstructuredGrid *grid,
  vtkIdType            Nnodes,
  bool                 overwrite
)
{
  EGVTKOBJECT_CREATENODEFIELD("node_status", vtkIntArray,    overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_layer",  vtkIntArray,    overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_index",  vtkLongArray_t, overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_meshdensity",  vtkDoubleArray, overwrite);
};

void EgVtkObject::allocateGrid
(
  vtkUnstructuredGrid *grid,
  vtkIdType            Ncells,
  vtkIdType            Nnodes
)
{
  EG_VTKSP(vtkPoints,points);
  points->SetNumberOfPoints(Nnodes);
  grid->SetPoints(points);
  grid->Allocate(Ncells,max(vtkIdType(1),Ncells/10));
  createBasicFields(grid, Ncells, Nnodes);
};

vec3_t EgVtkObject::cellCentre(vtkUnstructuredGrid *grid, vtkIdType id_cell)
{
  vec3_t x,xc(0,0,0);
  vtkIdType *pts, N_pts;
  grid->GetCellPoints(id_cell, N_pts, pts);
  double f = 1.0/N_pts;
  for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
    grid->GetPoint(pts[i_pts], x.data());
    xc += f*x;
  };
  return xc;
};

void EgVtkObject::getRestCells(vtkUnstructuredGrid      *grid, 
                               const QVector<vtkIdType> &cells,
                               QVector<vtkIdType>       &rest_cells)
{
  QVector<bool> is_in_cells(grid->GetNumberOfCells(), false);
  foreach (vtkIdType id_cell, cells) {
    is_in_cells[id_cell] = true;
  };
  rest_cells.resize(grid->GetNumberOfCells() - cells.size());
  int i_rest_cells = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (!is_in_cells[id_cell]) {
      rest_cells[i_rest_cells] = id_cell;
      ++i_rest_cells;
    };
  };
};

void EgVtkObject::makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst)
{
  allocateGrid(dst, src->GetNumberOfCells(), src->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    src->GetPoints()->GetPoint(id_node, x.data());
    dst->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(src, id_node, dst, id_node);
  };
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = src->GetCellType(id_cell);
    src->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType id_new_cell = dst->InsertNextCell(type_cell, N_pts, pts);
    copyCellData(src, id_cell, dst, id_new_cell);
  };
};

void EgVtkObject::makeCopyNoAlloc(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst)
{
  for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    src->GetPoints()->GetPoint(id_node, x.data());
    dst->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(src, id_node, dst, id_node);
  };
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = src->GetCellType(id_cell);
    src->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType id_new_cell = dst->InsertNextCell(type_cell, N_pts, pts);
    copyCellData(src, id_cell, dst, id_new_cell);
  };
};

int EgVtkObject::findVolumeCell
(
  vtkUnstructuredGrid      *grid,
  vtkIdType                id_surf,
  const QVector<int>       _nodes,      
  const QVector<vtkIdType> cells,      
  const QVector<int>       _cells,      
  QVector<QSet<int> >       &n2c
)
{
  vtkIdType N_pts, *pts;
  if (_cells.size()) N_pts = N_pts; // dummy statement to get rid of compiler warning ...
  grid->GetCellPoints(id_surf, N_pts, pts);
  QVector<QSet<int> > inters(N_pts-1);
  setIntersection(n2c[_nodes[pts[0]]], n2c[_nodes[pts[1]]], inters[0]);
  int i_pts = 2;
  while (i_pts < N_pts) {
    setIntersection(inters[i_pts-2], n2c[_nodes[pts[i_pts]]], inters[i_pts-1]);
    ++i_pts;
  };
  if (inters[N_pts-2].size() == 0) {
    return -1;
  } else if (inters[N_pts-2].size() > 2) {
    EG_BUG;
  };
  vtkIdType id_vol = -1;
  foreach (int i_cells, inters[N_pts-2]) {
    if (cells[i_cells] != id_surf) {
      id_vol = cells[i_cells];
    };
  };
  return id_vol;
};

void EgVtkObject::setBoundaryCodes(const QSet<int> &bcs)
{
  boundary_codes.clear();
  int bc;
  foreach(bc, bcs) {
    boundary_codes.insert(bc);
  };
};

void EgVtkObject::createIndices(vtkUnstructuredGrid *grid)
{
  if (!grid->GetCellData()->GetArray("cell_index")) {
    EG_VTKSP(vtkLongArray_t, var);
    var->SetName("cell_index");
    var->SetNumberOfValues(grid->GetNumberOfCells());
    grid->GetCellData()->AddArray(var);
  } else {
    EG_VTKDCC(vtkLongArray_t, var, grid, "cell_index");
    var->SetNumberOfValues(grid->GetNumberOfCells());
  };
  EG_VTKDCC(vtkLongArray_t, cell_index, grid, "cell_index");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cell_index->SetValue(id_cell, id_cell);
  };
  
  if (!grid->GetCellData()->GetArray("vtk_type")) {
    EG_VTKSP(vtkIntArray, var);
    var->SetName("vtk_type");
    var->SetNumberOfValues(grid->GetNumberOfCells());
    grid->GetCellData()->AddArray(var);
  } else {
    EG_VTKDCC(vtkIntArray, var, grid, "vtk_type");
    var->SetNumberOfValues(grid->GetNumberOfCells());
  };
  EG_VTKDCC(vtkIntArray, vtk_type, grid, "vtk_type");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtk_type->SetValue(id_cell, grid->GetCellType(id_cell));
  };
  
  if (!grid->GetCellData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, var);
    var->SetName("node_index");
    var->SetNumberOfValues(grid->GetNumberOfPoints());
    grid->GetPointData()->AddArray(var);
  } else {
    EG_VTKDCC(vtkLongArray_t, var, grid, "node_index");
    var->SetNumberOfValues(grid->GetNumberOfPoints());
  };
  EG_VTKDCN(vtkLongArray_t, node_index, grid, "node_index");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    node_index->SetValue(id_node, id_node);
  };
};

BoundaryCondition EgVtkObject::getBC(int bc)
{
  return GuiMainWindow::pointer()->getBC(bc);
};

int EgVtkObject::getSet(QString group, QString key, int value, int& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "int/" + key;
  qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toInt();
  qset->endGroup();
  return(variable);
}

double EgVtkObject::getSet(QString group, QString key, double value, double& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "double/" + key;
  qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toDouble();
  qset->endGroup();
  return(variable);
}

bool EgVtkObject::getSet(QString group, QString key, bool value, bool& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "bool/" + key;
  qset->beginGroup(group);
  Qt::CheckState state = (Qt::CheckState) ( value ? 2 : 0 );
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,state);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toBool();
  qset->endGroup();
  return(variable);
}

///////////////////////////////////////////
int cout_grid(ostream &stream, vtkUnstructuredGrid *grid, bool npoints, bool ncells, bool points, bool cells)
{
  stream<<"============="<<endl;
  if(npoints) stream << "grid->GetNumberOfPoints()=" << grid->GetNumberOfPoints() << endl;
  if(ncells) stream << "grid->GetNumberOfCells()=" << grid->GetNumberOfCells() << endl;
  if(points) {
    for (vtkIdType i = 0; i < grid->GetNumberOfPoints(); ++i) {
      vec3_t x;
      grid->GetPoint(i, x.data());
      stream << "Vertex " << i << " = " << x << endl;
    };
  }
  if(cells) {
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    for (vtkIdType i = 0; i < grid->GetNumberOfCells(); ++i) {
      vtkCell *C = (vtkCell *) vtkCell::New();
      C=grid->GetCell(i);
      vtkIdType npts=C->GetNumberOfPoints();
      vtkIdType* pts;
      grid->GetCellPoints(i, npts, pts);
      stream << "Cell " << i << " = ";
      for(int j=0;j<npts;j++) stream << pts[j] << " ";
      stream << "boundary_code=" << cell_code->GetValue(i);
      stream << endl;
    };
  }
  stream<<"============="<<endl;
  return 0;
}

///////////////////////////////////////////
int addPoint(vtkUnstructuredGrid* a_grid,vtkIdType index,vec3_t x)
{
  a_grid->GetPoints()->SetPoint(index, x.data());
  return(0);
}

int addCell(vtkUnstructuredGrid* a_grid, vtkIdType A, vtkIdType B, vtkIdType C, int bc)
{
  vtkIdType npts=3;
  vtkIdType pts[3];
  pts[0]=A;
  pts[1]=B;
  pts[2]=C;
  vtkIdType newCellId = a_grid->InsertNextCell(VTK_TRIANGLE,npts,pts);
  EG_VTKDCC(vtkIntArray, cell_code, a_grid, "cell_code");
  cell_code->SetValue(newCellId, bc);
  return(0);
}

///////////////////////////////////////////

int getShortestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_id_cell, N_pts, pts);
  vec3_t* x=new vec3_t[N_pts];
  for(int i=0;i<N_pts;i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_minlen=0;
  double minlen=(x[1]-x[0]).abs();
  for(int i=1;i<N_pts;i++)
  {
    double len=(x[(i+1)%N_pts]-x[i]).abs();
    if(len<minlen){
      minlen=len;
      id_minlen=i;
    }
  }
  delete x;
  return(id_minlen);
}

int getLongestSide(vtkIdType a_id_cell,vtkUnstructuredGrid* a_grid)
{
  vtkIdType N_pts, *pts;
  a_grid->GetCellPoints(a_id_cell, N_pts, pts);
  vec3_t* x=new vec3_t[N_pts];
  for(int i=0;i<N_pts;i++) a_grid->GetPoints()->GetPoint(pts[i], x[i].data());
  int id_maxlen=0;
  double maxlen=(x[1]-x[0]).abs();
  cout<<"maxlen="<<maxlen<<endl;
  for(int i=1;i<N_pts;i++)
  {
    double len=(x[(i+1)%N_pts]-x[i]).abs();
    cout<<"len["<<i<<"]="<<len<<endl;
    if(len>maxlen){
      maxlen=len;
      id_maxlen=i;
    }
  }
  delete x;
  return(id_maxlen);
}
///////////////////////////////////////////

QSet <int> complementary_bcs(QSet <int> &bcs, vtkUnstructuredGrid *a_grid, QVector <vtkIdType> &a_cells)
{
  QSet <int> bcs_complement;
  EG_VTKDCC(vtkIntArray, bc, a_grid, "cell_code");
  foreach (vtkIdType id_cell, a_cells) {
    int code=bc->GetValue(id_cell);
    if (!bcs.contains(code)) {
      bcs_complement.insert(code);
    };
  }
  return(bcs_complement);
}

///////////////////////////////////////////

QString cell2str(vtkIdType id_cell,vtkUnstructuredGrid* grid)
{
  QString tmp;
  tmp.setNum(id_cell);
  QString txt = "id_cell=" + tmp;
  
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(id_cell, N_pts, pts);
  
  txt += " [";
  for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
    tmp.setNum(pts[i_pts]);
    txt += tmp;
    if (i_pts < N_pts-1) {
      txt += ",";
    };
  };
  txt += "]";
  return(txt);
}



///////////////////////////////////////////

Qt::CheckState int2CheckState(int a)
{
  if(a==0) return(Qt::Unchecked);
//   if(a==1) return(Qt::PartiallyChecked);
  else return(Qt::Checked);
}

// ///////////////////////////////////////////
// /* Here is how we we get QTextStreams that look like iostreams */
// Qcin=QTextStream(stdin, QIODevice::ReadOnly);
// Qcout=QTextStream(stdout, QIODevice::WriteOnly);
// Qcerr=QTextStream(stderr, QIODevice::WriteOnly);
// ///////////////////////////////////////////

pair<vtkIdType,vtkIdType> OrderedPair(vtkIdType a, vtkIdType b)
{
  vtkIdType x=min(a,b);
  vtkIdType y=max(a,b);
  return(pair<vtkIdType,vtkIdType>(x,y));
}
