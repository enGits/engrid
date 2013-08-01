// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
#include <vtkCharArray.h>
#include <vtkCellArray.h>

int EgVtkObject::DebugLevel;

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
      }
    }
  }
  for (int i_node = 0; i_node < nodes.count(); ++i_node) {
    node_normals[i_node].normalise();
    //cout << node_normals[i_node] << endl;
  }
}

void EgVtkObject::createNodeMapping
(
  QVector<vtkIdType>  &nodes,
  QVector<int>        &_nodes,
  vtkUnstructuredGrid *grid
)
{
  _nodes.fill(-1,grid->GetNumberOfPoints());
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    _nodes[nodes[i_nodes]] = i_nodes;
  }
}

void EgVtkObject::createCellMapping
(
  QVector<vtkIdType>  &cells,
  QVector<int>        &_cells,
  vtkUnstructuredGrid *grid
)
{
  _cells.fill(-1,grid->GetNumberOfCells());
  for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
    _cells[cells[i_cells]] = i_cells;
  }
}

void EgVtkObject::createNodeToBcMapping
(
  QVector<QSet<int> >  &bcs,
  vtkUnstructuredGrid  *grid
)
{
  EG_BUG;
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
        }
      }
    }
  }
}

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
    QList<vtkIdType> pts;
    getPointsOfCell(grid, cells[i_cells], pts);
    foreach (vtkIdType id_node, pts) {
      n2c[_nodes[id_node]].insert(i_cells);
    }
  }
}

void EgVtkObject::createNodeToCell
(
  QVector<vtkIdType>     &cells,
  QVector<vtkIdType>     &nodes,
  QVector<int>           &_nodes,
  QVector<QVector<int> > &n2c,
  vtkUnstructuredGrid    *grid
)
{
  n2c.fill(QVector<int>(), nodes.size());
  QVector<int> count(nodes.size(),0);
  for (vtkIdType i_cells = 0; i_cells < cells.size(); ++i_cells) {
    QList<vtkIdType> pts;
    getPointsOfCell(grid, cells[i_cells], pts);
    foreach (vtkIdType id_node, pts) {
      ++count[_nodes[id_node]];
    }
  }
  for (int i = 0; i < nodes.size(); ++i) {
    n2c[i].resize(count[i]);
    count[i] = 0;
  }
  for (vtkIdType i_cells = 0; i_cells < cells.size(); ++i_cells) {
    QList<vtkIdType> pts;
    getPointsOfCell(grid, cells[i_cells], pts);
    foreach (vtkIdType id_node, pts) {
      int i_nodes = _nodes[id_node];
      n2c[i_nodes][count[i_nodes]] = i_cells;
      ++count[i_nodes];
    }
  }
}

void EgVtkObject::addToN2N(QVector<QSet<int> > &n2n, int n1, int n2)
{
  n2n[n1].insert(n2);
  n2n[n2].insert(n1);
}

void EgVtkObject::createNodeToNode(QVector<vtkIdType> &cells, QVector<vtkIdType> &nodes, QVector<int> &_nodes, QVector<QSet<int> > &n2n, vtkUnstructuredGrid *grid)
{
  n2n.fill(QSet<int>(), nodes.size());
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, grid)) {
      vtkIdType num_pts, *pts;
      grid->GetCellPoints(id_cell, num_pts, pts);
      QVector<int> n(num_pts);
      for (int i = 0; i < num_pts; ++i) {
        n[i] = _nodes[pts[i]];
      }
      for (int i = 0; i < num_pts - 1; ++i) {
        addToN2N(n2n, n[i], n[i+1]);
      }
      addToN2N(n2n, n.last(), n.first());
    } else {
      vtkIdType num_faces = 0;
      vtkIdType type_cell = grid->GetCellType(id_cell);

      if        (type_cell == VTK_TETRA) {
        num_faces = 4;
      } else if (type_cell == VTK_WEDGE) {
        num_faces = 5;
      } else if (type_cell == VTK_PYRAMID) {
        num_faces = 5;
      } else if (type_cell == VTK_HEXAHEDRON) {
        num_faces = 6;
      } else if (type_cell == VTK_POLYHEDRON) {
        vtkIdType *pts;
        grid->GetFaceStream(id_cell, num_faces, pts);
      }

      for (int i = 0; i < num_faces; ++i) {
        QVector<vtkIdType> face;
        getFaceOfCell(grid, id_cell, i, face);
        QVector<int> n(face.size());
        for (int j = 0; j < face.size(); ++j) {
          n[j] = _nodes[face[j]];
        }
        for (int j = 0; j < face.size() - 1; ++j) {
          addToN2N(n2n, n[j], n[j+1]);
        }
        addToN2N(n2n, n.last(), n.first());
      }
    }



    /*
    vtkIdType *pts;
    vtkIdType  Npts;
    grid->GetCellPoints(id_cell, Npts, pts);
    vector<int> n(Npts);
    for (int i = 0; i < Npts; ++i) {
      n[i] = _nodes[pts[i]];
    }
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_TRIANGLE) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[2], n[0]);
    } else if (type_cell == VTK_QUAD) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[2], n[3]);
      addToN2N(n2n, n[3], n[0]);
    } else if (type_cell == VTK_POLYGON) {
      EG_BUG;
    } else if (type_cell == VTK_TETRA) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[2]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[3]);
      addToN2N(n2n, n[2], n[3]);
    } else if (type_cell == VTK_PYRAMID) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[0], n[4]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[4]);
      addToN2N(n2n, n[2], n[3]);
      addToN2N(n2n, n[2], n[4]);
      addToN2N(n2n, n[3], n[4]);
    } else if (type_cell == VTK_WEDGE) {
      addToN2N(n2n, n[0], n[1]);
      addToN2N(n2n, n[0], n[2]);
      addToN2N(n2n, n[0], n[3]);
      addToN2N(n2n, n[1], n[2]);
      addToN2N(n2n, n[1], n[4]);
      addToN2N(n2n, n[2], n[5]);
      addToN2N(n2n, n[3], n[4]);
      addToN2N(n2n, n[3], n[5]);
      addToN2N(n2n, n[4], n[5]);
    } else if (type_cell == VTK_HEXAHEDRON) {
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
    } else if (type_cell == VTK_POLYHEDRON) {
      EG_BUG;
    }
    */
  }
}

void EgVtkObject::createNodeToNode
(
  QVector<vtkIdType>     &cells,
  QVector<vtkIdType>     &nodes,
  QVector<int>           &_nodes,
  QVector<QVector<int> > &n2n,
  vtkUnstructuredGrid    *grid
)
{
  QVector<QSet<int> > n2n_set;
  createNodeToNode(cells, nodes, _nodes, n2n_set, grid);
  n2n.resize(n2n_set.size());
  for (int i = 0; i < n2n.size(); ++i) {
    n2n[i].resize(n2n_set[i].size());
    qCopy(n2n_set[i].begin(), n2n_set[i].end(), n2n[i].begin());
  }
}

void EgVtkObject::getAllCells
(
  QVector<vtkIdType>  &cells,
  vtkUnstructuredGrid *grid
)
{
  int N = 0;
  cells.resize(grid->GetNumberOfCells());
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cells[N] = id_cell;
    ++N;
  }
}

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
    }
  }
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (grid->GetCellType(id_cell) == type) {
      cells[N] = id_cell;
      ++N;
    }
  }
}


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
    }
  }
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, grid)) {
      cells[N] = id_cell;
      ++N;
    }
  }
}

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
    }
  }
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      cells[N] = id_cell;
      ++N;
    }
  }
}

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
      }
    }
  }
  cells.resize(N);
  N = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      if (bcs.contains(cell_code->GetValue(id_cell))) {
        cells[N] = id_cell;
        ++N;
      }
    }
  }
}

void EgVtkObject::addToC2C(vtkIdType id_cell, QVector<int> &_cells, QVector<QVector<int> > &c2c, int j, vtkIdList *nds, vtkIdList *cls, vtkUnstructuredGrid *grid)
{
  c2c[_cells[id_cell]][j] = -1;
  grid->GetCellNeighbors(id_cell, nds, cls);
  if (isSurface(id_cell, grid)) {
    for (int i = 0; i < cls->GetNumberOfIds(); ++i) {
      if (cls->GetId(i) != id_cell) {
        if (_cells[cls->GetId(i)] != -1) {
          if (isSurface(cls->GetId(i), grid)) {
            c2c[_cells[id_cell]][j] = _cells[cls->GetId(i)];
          }
        }
      }
    }
  } else {
    for (int i = 0; i < cls->GetNumberOfIds(); ++i) {
      if (cls->GetId(i) != id_cell) {
        if (_cells[cls->GetId(i)] != -1) {
          if (isVolume(cls->GetId(i), grid) || c2c[_cells[id_cell]][j] == -1) {
            c2c[_cells[id_cell]][j] = _cells[cls->GetId(i)];
          }
        }
      }
    }
  }
}


void EgVtkObject::createCellToCell(QVector<vtkIdType> &cells, QVector<QVector<int> > &c2c, vtkUnstructuredGrid *grid)
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
    //vtkIdType *pts;
    //vtkIdType  Npts;
    //grid->GetCellPoints(id_cell, Npts, pts);
    QList<vtkIdType> pts;
    getPointsOfCell(grid, id_cell, pts);
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
    } else if (grid->GetCellType(id_cell) == VTK_POLYGON) {
      c2c[i].resize(pts.size());
      pts.append(pts[0]);
      for (int i = 0; i < pts.size() - 1; ++i) {
        nds->Reset();
        nds->InsertNextId(pts[i]);
        nds->InsertNextId(pts[i+1]);
        addToC2C(id_cell, _cells, c2c, i, nds, cls, grid);
      }
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
    } else if (grid->GetCellType(id_cell) == VTK_POLYHEDRON) {
      EG_VTKSP(vtkIdList, stream);
      grid->GetFaceStream(id_cell, stream);
      int num_faces = stream->GetId(0);
      c2c[i].resize(num_faces);
      int id = 1;
      for (int i_face = 0; i_face < num_faces; ++i_face) {
        nds->Reset();
        int num_pts = stream->GetId(id++);
        for (int i_pt = 0; i_pt < num_pts; ++i_pt) {
          nds->InsertNextId(stream->GetId(id++));
        }
        addToC2C(id_cell, _cells, c2c, i_face, nds, cls, grid);
      }
    }
  }
}

bool EgVtkObject::isVolume(vtkIdType id_cell, vtkUnstructuredGrid *grid)
{
  bool isVol = false;
  if      (grid->GetCellType(id_cell) == VTK_TETRA)      isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_PYRAMID)    isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_WEDGE)      isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_HEXAHEDRON) isVol = true;
  else if (grid->GetCellType(id_cell) == VTK_POLYHEDRON) isVol = true;
  return isVol;
}

bool EgVtkObject::isSurface(vtkIdType id_cell, vtkUnstructuredGrid *grid)
{
  bool isSurf = false;
  if      (grid->GetCellType(id_cell) == VTK_TRIANGLE) isSurf = true;
  else if (grid->GetCellType(id_cell) == VTK_QUAD)     isSurf = true;
  else if (grid->GetCellType(id_cell) == VTK_POLYGON)  isSurf = true;
  return isSurf;
}

void EgVtkObject::UpdateCellIndex(vtkUnstructuredGrid *grid)
{
  if (!grid->GetCellData()->GetArray("cell_index")) {
    EG_VTKSP(vtkLongArray_t, cell_index);
    cell_index->SetName("cell_index");
    cell_index->SetNumberOfValues(grid->GetNumberOfCells());
    grid->GetCellData()->AddArray(cell_index);
  }
  EG_VTKDCC(vtkLongArray_t, cell_index, grid, "cell_index");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cell_index->SetValue(id_cell, id_cell);
  }
}

void EgVtkObject::UpdateNodeIndex(vtkUnstructuredGrid *grid)
{
  if (!grid->GetPointData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, node_index);
    node_index->SetName("node_index");
    node_index->SetNumberOfValues(grid->GetNumberOfPoints());
    grid->GetPointData()->AddArray(node_index);
  }
  EG_VTKDCN(vtkLongArray_t, node_index, grid, "node_index");
  for (vtkIdType pointId = 0; pointId < grid->GetNumberOfPoints(); ++pointId) {
    node_index->SetValue(pointId, pointId);
  }
}

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
  }
  if (!pdata->GetPointData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, node_index);
    node_index->SetName("node_index");
    //node_index->SetNumberOfValues(nodes.size());
    pdata->GetPointData()->AddArray(node_index);
  }
  EG_VTKDCC(vtkLongArray_t, pd_cell_index, pdata, "cell_index");
  EG_VTKDCN(vtkLongArray_t, pd_node_index, pdata, "node_index");
  pd_cell_index->SetNumberOfValues(cells.size());
  pd_node_index->SetNumberOfValues(nodes.size());
  for (int i_cell = 0; i_cell < cells.size(); ++i_cell) {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType cellType = grid->GetCellType(id_cell);
    if ((cellType != VTK_TRIANGLE) && (cellType != VTK_QUAD)) {
      EG_ERR_RETURN("unsupported cell type for this operation");
    }
    vtkIdType Npts, *pts;
    grid->GetCellPoints(id_cell, Npts, pts);
    vtkIdType *new_pts = new vtkIdType[Npts];
    for (int i = 0; i < Npts; ++i) {
      new_pts[i] = _nodes[pts[i]];
    }
    vtkIdType newCellId = pdata->InsertNextCell(cellType, Npts, new_pts);
    pd_cell_index->SetValue(newCellId, id_cell);
    delete [] new_pts;
  }
  for (int i_node = 0; i_node < nodes.size(); ++i_node) {
    vec3_t x;
    grid->GetPoints()->GetPoint(nodes[i_node], x.data());
    pdata->GetPoints()->SetPoint(i_node, x.data());
    pd_node_index->SetValue(i_node, nodes[i_node]);
  }
}

#define EGVTKOBJECT_COPYCELLDATA(FIELD,TYPE) \
{ \
if (old_grid->GetCellData()->GetArray(FIELD)) { \
EG_VTKDCC(TYPE, var1, old_grid, FIELD); \
EG_VTKDCC(TYPE, var2, new_grid, FIELD); \
var2->SetValue(newId, var1->GetValue(oldId)); \
} \
}

void EgVtkObject::copyCellData
(
  vtkUnstructuredGrid *old_grid, 
  vtkIdType            oldId, 
  vtkUnstructuredGrid *new_grid,
  vtkIdType            newId
)
{
  EGVTKOBJECT_COPYCELLDATA("vtk_type",    vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_code",   vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_orgdir", vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_curdir", vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_voldir", vtkIntArray);
  EGVTKOBJECT_COPYCELLDATA("cell_index",  vtkLongArray_t);
}

#define EGVTKOBJECT_COPYNODEDATA(FIELD,TYPE) \
{ \
if (old_grid->GetPointData()->GetArray(FIELD)) { \
EG_VTKDCN(TYPE, var1, old_grid, FIELD); \
EG_VTKDCN(TYPE, var2, new_grid, FIELD); \
var2->SetValue(newId, var1->GetValue(oldId)); \
} \
}

void EgVtkObject::copyNodeData
(
  vtkUnstructuredGrid *old_grid, 
  vtkIdType            oldId, 
  vtkUnstructuredGrid *new_grid,
  vtkIdType            newId
)
{
  EGVTKOBJECT_COPYNODEDATA("node_status", vtkIntArray);
  EGVTKOBJECT_COPYNODEDATA("node_layer", vtkIntArray);
  EGVTKOBJECT_COPYNODEDATA("node_index", vtkLongArray_t);
  EGVTKOBJECT_COPYNODEDATA("node_specified_density", vtkIntArray);
  EGVTKOBJECT_COPYNODEDATA("node_meshdensity_desired",  vtkDoubleArray);
  //EGVTKOBJECT_COPYNODEDATA("node_meshdensity_current",  vtkDoubleArray);
  EGVTKOBJECT_COPYNODEDATA("node_type",  vtkCharArray);
  EGVTKOBJECT_COPYNODEDATA("node_pindex", vtkLongArray_t);
}

#define EGVTKOBJECT_CREATECELLFIELD(FIELD,TYPE,OW) \
if (!grid->GetCellData()->GetArray(FIELD)) { \
  EG_VTKSP(TYPE, var); \
  var->SetName(FIELD); \
  var->SetNumberOfValues(Ncells); \
  grid->GetCellData()->AddArray(var); \
  for (int i = 0; i < grid->GetNumberOfCells(); ++i) { \
    var->SetValue(i,0); \
  } \
} else if (OW) { \
  EG_VTKDCC(TYPE, var, grid, FIELD); \
  var->SetNumberOfValues(Ncells); \
  for (int i = 0; i < grid->GetNumberOfCells(); ++i) { \
    var->SetValue(i,0); \
  } \
}

#define EGVTKOBJECT_CREATENODEFIELD(FIELD,TYPE,OW) \
if (!grid->GetPointData()->GetArray(FIELD)) { \
  EG_VTKSP(TYPE, var); \
  var->SetName(FIELD); \
  var->SetNumberOfValues(Nnodes); \
  grid->GetPointData()->AddArray(var); \
  for (int i = 0; i < grid->GetNumberOfPoints(); ++i) { \
    var->SetValue(i,0); \
  } \
} else if (OW) { \
  EG_VTKDCN(TYPE, var, grid, FIELD); \
  var->SetNumberOfValues(Nnodes); \
  for (int i = 0; i < grid->GetNumberOfPoints(); ++i) { \
    var->SetValue(i,0); \
  } \
}

void EgVtkObject::createBasicFields(vtkUnstructuredGrid *grid, vtkIdType Ncells, vtkIdType Nnodes, bool overwrite)
{
  createBasicNodeFields(grid, Nnodes, overwrite);
  createBasicCellFields(grid, Ncells, overwrite);
}

void EgVtkObject::createBasicCellFields(vtkUnstructuredGrid *grid, vtkIdType Ncells, bool overwrite)
{
  EGVTKOBJECT_CREATECELLFIELD("vtk_type" ,         vtkIntArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_code",         vtkIntArray, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_index",        vtkLongArray_t, overwrite);
  EGVTKOBJECT_CREATECELLFIELD("cell_orgdir",       vtkIntArray, overwrite); // original orientation
  EGVTKOBJECT_CREATECELLFIELD("cell_curdir",       vtkIntArray, overwrite); // current orientation
  EGVTKOBJECT_CREATECELLFIELD("cell_voldir",       vtkIntArray, overwrite); // volume orientation -- only valid for a single (i.e. the current) volume
  EGVTKOBJECT_CREATECELLFIELD("cell_mesh_quality", vtkDoubleArray, overwrite); // generic field to store different quality measures
}

void EgVtkObject::createBasicNodeFields(vtkUnstructuredGrid *grid, vtkIdType Nnodes, bool overwrite)
{
  EGVTKOBJECT_CREATENODEFIELD("node_status",               vtkIntArray,    overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_layer",                vtkIntArray,    overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_index",                vtkLongArray_t, overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_specified_density",    vtkIntArray,    overwrite); //density index from table
  EGVTKOBJECT_CREATENODEFIELD("node_meshdensity_desired",  vtkDoubleArray, overwrite); //what we want
  EGVTKOBJECT_CREATENODEFIELD("node_type",                 vtkCharArray,   overwrite); //node type
  EGVTKOBJECT_CREATENODEFIELD("node_type_counter",         vtkIntArray,    overwrite); // counter field to delay node type demotion
  EGVTKOBJECT_CREATENODEFIELD("node_pindex",               vtkLongArray_t, overwrite);
  EGVTKOBJECT_CREATENODEFIELD("node_mesh_quality",         vtkDoubleArray, overwrite); // generic field to store different quality measures
}

void EgVtkObject::allocateGrid(vtkUnstructuredGrid *grid, vtkIdType Ncells, vtkIdType Nnodes, bool create_fields)
{
  EG_VTKSP(vtkPoints, points);
  points->SetNumberOfPoints(Nnodes);
  grid->SetPoints(points);
  grid->Allocate(0, 0);
  grid->Squeeze();
  grid->Reset();
  grid->Allocate(Ncells,max(vtkIdType(1),Ncells/10));
  if (create_fields) {
    createBasicFields(grid, Ncells, Nnodes, true);
  }
}

vec3_t EgVtkObject::cellCentre(vtkUnstructuredGrid *grid, vtkIdType id_cell)
{
  vec3_t x,xc(0,0,0);
  vtkIdType *pts, N_pts;
  grid->GetCellPoints(id_cell, N_pts, pts);
  double f = 1.0/N_pts;
  for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
    grid->GetPoint(pts[i_pts], x.data());
    xc += f*x;
  }
  return xc;
}

double EgVtkObject::faceAngle(vtkUnstructuredGrid *grid, vtkIdType id_face1, vtkIdType id_face2)
{
  QList<vtkIdType> edge_nodes;
  sharedNodesOfCells(grid, id_face1, id_face2, edge_nodes);
  if (edge_nodes.size() != 2) {
    EG_BUG;
  }
  vec3_t x1, x2;
  grid->GetPoint(edge_nodes[0], x1.data());
  grid->GetPoint(edge_nodes[1], x2.data());
  vec3_t xf1 = cellCentre(grid, id_face1);
  vec3_t xf2 = cellCentre(grid, id_face2);
  vec3_t xe  = 0.5*(x1 + x2);
  vec3_t xfe = 0.5*(xf1 + xf2);
  vec3_t n1  = cellNormal(grid, id_face1);
  vec3_t n2  = cellNormal(grid, id_face2);
  vec3_t ve  = xe - xfe;
  double alpha = angle(n1, n2);
  if ((n1 + n2)*ve < 0) {
    return -alpha;
  }
  return alpha;
}

void EgVtkObject::getRestCells(vtkUnstructuredGrid      *grid, 
                               const QVector<vtkIdType> &cells,
                               QVector<vtkIdType>       &rest_cells)
{
  QVector<bool> is_in_cells(grid->GetNumberOfCells(), false);
  foreach (vtkIdType id_cell, cells) {
    is_in_cells[id_cell] = true;
  }
  rest_cells.resize(grid->GetNumberOfCells() - cells.size());
  int i_rest_cells = 0;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (!is_in_cells[id_cell]) {
      rest_cells[i_rest_cells] = id_cell;
      ++i_rest_cells;
    }
  }
}

void EgVtkObject::makeCopy(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, bool copy_data)
{
  allocateGrid(dst, src->GetNumberOfCells(), src->GetNumberOfPoints(), copy_data);
  for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    src->GetPoints()->GetPoint(id_node, x.data());
    dst->GetPoints()->SetPoint(id_node, x.data());
    if (copy_data) {
      copyNodeData(src, id_node, dst, id_node);
    }
  }
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
    vtkIdType id_new_cell = copyCell(src, id_cell, dst);
    /*
    if (id_cell == 2995) {
      vtkUnstructuredGrid *G = dst;
      EG_VTKSP(vtkIdList, stream);
      vtkIdType type_cell = G->GetCellType(id_cell);
      if (type_cell == VTK_POLYHEDRON) {
        G->GetFaceStream(id_cell, stream);
      } else {
        G->GetCellPoints(id_cell, stream);
      }
      QList<int> ids;
      vtkIdList2QContainer(stream, ids);
      vtkIdType *pointer = stream->GetPointer(18);
      cerr << "break" << endl;
    }
    */
    if (copy_data) {
      copyCellData(src, id_cell, dst, id_new_cell);
    }
  }
  /*
  if (!dbg) return;

  {
    vtkIdType id_cell = 2995;
    vtkUnstructuredGrid *G = src;
    EG_VTKSP(vtkIdList, stream);
    vtkIdType type_cell = G->GetCellType(id_cell);
    if (type_cell == VTK_POLYHEDRON) {
      G->GetFaceStream(id_cell, stream);
    } else {
      G->GetCellPoints(id_cell, stream);
    }
    QList<int> ids;
    vtkIdList2QContainer(stream, ids);
    vtkIdType *pointer = stream->GetPointer(18);
    cerr << "break" << endl;
  }
  {
    vtkIdType id_cell = 2995;
    vtkUnstructuredGrid *G = dst;
    EG_VTKSP(vtkIdList, stream);
    vtkIdType type_cell = G->GetCellType(id_cell);
    if (type_cell == VTK_POLYHEDRON) {
      G->GetFaceStream(id_cell, stream);
    } else {
      G->GetCellPoints(id_cell, stream);
    }
    QList<int> ids;
    vtkIdList2QContainer(stream, ids);
    vtkIdType *pointer = stream->GetPointer(18);
    cerr << "break" << endl;
  }
  */
}

void EgVtkObject::makeCopyNoAlloc(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst)
{
  for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    src->GetPoints()->GetPoint(id_node, x.data());
    dst->GetPoints()->SetPoint(id_node, x.data());
    copyNodeData(src, id_node, dst, id_node);
  }
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
    vtkIdType num, *stream;
    vtkIdType type_cell = src->GetCellType(id_cell);
    if (type_cell == VTK_POLYHEDRON) {
      src->GetFaceStream(id_cell, num, stream);
    } else {
      src->GetCellPoints(id_cell, num, stream);
    }
    vtkIdType id_new_cell = dst->InsertNextCell(type_cell, num, stream);
    copyCellData(src, id_cell, dst, id_new_cell);
  }
}

void EgVtkObject::reorientateFace(vtkUnstructuredGrid *grid, vtkIdType id_face)
{
  EG_VTKDCC(vtkIntArray, cell_curdir, grid, "cell_curdir");
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(id_face, N_pts, pts);
  QVector<vtkIdType> new_pts(N_pts);
  for (int i = 0; i < N_pts; ++i) {
    new_pts[i] = pts[N_pts - i - 1];
  }
  if (cell_curdir->GetValue(id_face) == 0) {
    cell_curdir->SetValue(id_face, 1);
  } else {
    cell_curdir->SetValue(id_face, 0);
  }
  grid->ReplaceCell(id_face, N_pts, new_pts.data());
}

void EgVtkObject::resetOrientation(vtkUnstructuredGrid *grid)
{
  EG_VTKDCC(vtkIntArray, cell_orgdir, grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_curdir, grid, "cell_curdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, grid, "cell_voldir");
  QVector<vtkIdType> faces;
  getAllSurfaceCells(faces, grid);
  foreach (vtkIdType id_face, faces) {
    if (cell_curdir->GetValue(id_face) != cell_orgdir->GetValue(id_face)) {
      reorientateFace(grid, id_face);
      cell_curdir->SetValue(id_face, cell_orgdir->GetValue(id_face));
    }
    cell_voldir->SetValue(id_face, 0);
  }
}

vtkIdType EgVtkObject::findVolumeCell(vtkUnstructuredGrid *grid, vtkIdType id_surf, g2l_t _nodes, l2g_t cells, g2l_t _cells, l2l_t n2c)
{
  vtkIdType N_pts=0, *pts=NULL; //allways initialize variables, otherwise unexpected results will occur!
  grid->GetCellPoints(id_surf, N_pts, pts);
  QVector<QSet<int> > inters(N_pts-1);
  qcontIntersection(n2c[_nodes[pts[0]]], n2c[_nodes[pts[1]]], inters[0]);
  int i_pts = 2;
  while (i_pts < N_pts) {
    qcontIntersection(inters[i_pts-2], n2c[_nodes[pts[i_pts]]], inters[i_pts-1]);
    ++i_pts;
  }
  if (inters[N_pts-2].size() == 0) {
    return -1;
  } else if (inters[N_pts-2].size() > 2) {
    EG_BUG;
  }
  vtkIdType id_vol = -1;
  foreach (int i_cells, inters[N_pts-2]) {
    if (cells[i_cells] != id_surf) {
      id_vol = cells[i_cells];
    }
  }
  return id_vol;
}

void EgVtkObject::setBoundaryCodes(const QSet<int> &bcs)
{
  m_BoundaryCodes = bcs;
}

QSet<int> EgVtkObject::getBoundaryCodes()
{
  return m_BoundaryCodes;
}

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
  }
  EG_VTKDCC(vtkLongArray_t, cell_index, grid, "cell_index");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    cell_index->SetValue(id_cell, id_cell);
  }
  
  if (!grid->GetCellData()->GetArray("vtk_type")) {
    EG_VTKSP(vtkIntArray, var);
    var->SetName("vtk_type");
    var->SetNumberOfValues(grid->GetNumberOfCells());
    grid->GetCellData()->AddArray(var);
  } else {
    EG_VTKDCC(vtkIntArray, var, grid, "vtk_type");
    var->SetNumberOfValues(grid->GetNumberOfCells());
  }
  EG_VTKDCC(vtkIntArray, vtk_type, grid, "vtk_type");
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtk_type->SetValue(id_cell, grid->GetCellType(id_cell));
  }
  
  if (!grid->GetCellData()->GetArray("node_index")) {
    EG_VTKSP(vtkLongArray_t, var);
    var->SetName("node_index");
    var->SetNumberOfValues(grid->GetNumberOfPoints());
    grid->GetPointData()->AddArray(var);
  } else {
    EG_VTKDCC(vtkLongArray_t, var, grid, "node_index");
    var->SetNumberOfValues(grid->GetNumberOfPoints());
  }
  EG_VTKDCN(vtkLongArray_t, node_index, grid, "node_index");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    node_index->SetValue(id_node, id_node);
  }
}

BoundaryCondition EgVtkObject::getBC(int bc)
{
  return GuiMainWindow::pointer()->getBC(bc);
}

int EgVtkObject::getSet(QString group, QString key, int value, int& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "int/" + key;
  if(group!=QObject::tr("General")) qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toInt();
  if(group!=QObject::tr("General")) qset->endGroup();
  return(variable);
}

double EgVtkObject::getSet(QString group, QString key, double value, double& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "double/" + key;
  if(group!=QObject::tr("General")) qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toDouble();
  if(group!=QObject::tr("General")) qset->endGroup();
  return(variable);
}

bool EgVtkObject::getSet(QString group, QString key, bool value, bool& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key = "bool/" + key;
  if(group!=QObject::tr("General")) qset->beginGroup(group);
  Qt::CheckState state = (Qt::CheckState) ( value ? 2 : 0 );
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key,state);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key,variable)).toBool();
  if(group!=QObject::tr("General")) qset->endGroup();
  return(variable);
}

QString EgVtkObject::getSet(QString group, QString key, QString value, QString& variable)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key;
  typed_key = QObject::tr("QString/") + key;
  if (group != QObject::tr("General")) qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key, value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key)).toString();
  if (group != QObject::tr("General")) qset->endGroup();
  return(variable);
}

QString EgVtkObject::getSet(QString group, QString key, QString value, QString& variable, int type)
{
  QSettings *qset = GuiMainWindow::settings();
  QString typed_key;
  if (type == 0) {
    typed_key = QObject::tr("QString/") + key;
  }
  else if (type == 1) {
    typed_key = QObject::tr("Filename/") + key;
  }
  else {
    typed_key = QObject::tr("Directory/") + key;
  }
  if (group != QObject::tr("General")) qset->beginGroup(group);
  //if key=value pair not found in settings file, write it
  if (!qset->contains(typed_key)) qset->setValue(typed_key, value);
  //read key value from settings file and assign it to variable
  variable = (qset->value(typed_key)).toString();
  if (group != QObject::tr("General")) qset->endGroup();
  return(variable);
}

void EgVtkObject::writeGrid(vtkUnstructuredGrid *grid, QString name)
{
  QVector<vtkIdType> cells;
  getAllCells(cells, grid);
  name = GuiMainWindow::pointer()->getCwd() + "/" + name + ".vtu";
  writeCells(grid, cells, name);
}

void EgVtkObject::getAllNodeDataNames(QVector<QString> &field_names, vtkUnstructuredGrid *grid)
{
  int N = grid->GetPointData()->GetNumberOfArrays();
  field_names.resize(N);
  for (int i = 0; i < N; ++i) {
    field_names[i] = grid->GetPointData()->GetArrayName(i);
  }
}

void EgVtkObject::getAllCellDataNames(QVector<QString> &field_names, vtkUnstructuredGrid *grid)
{
  int N = grid->GetCellData()->GetNumberOfArrays();
  field_names.resize(N);
  for (int i = 0; i < N; ++i) {
    field_names[i] = grid->GetCellData()->GetArrayName(i);
  }
}

QString EgVtkObject::stripFromExtension(QString file_name)
{
  int i = file_name.size() - 1;
  while ((i > 0) && (file_name[i] != '.') && (file_name[i] != '/') && (file_name[i] != '\\')) {
    --i;
  }
  if (file_name[i] == '.') {
    return file_name.left(i);
  }
  return file_name;
}

QString EgVtkObject::getExtension(QString file_name)
{
  int i = file_name.size();
  while ((i > 0) && (file_name[i] != '.') && (file_name[i] != '/') && (file_name[i] != '\\')) {
    --i;
  }
  if (file_name[i] == '.') {
    return (file_name.right(file_name.size() - i - 1)).toLower();
  }
  return "";
}

///////////////////////////////////////////

void EgVtkObject::getFaceOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face, QVector<vtkIdType> &ids)
{
  vtkIdType type_cell = grid->GetCellType(id_cell);
  ids.clear();
  EG_VTKSP(vtkIdList, stream);
  QVector<vtkIdType> pts;
  if (type_cell == VTK_POLYHEDRON) {
    grid->GetFaceStream(id_cell, stream);
  } else {
    grid->GetCellPoints(id_cell, stream);
  }
  vtkIdList2QContainer(stream, pts);

  if (type_cell == VTK_TETRA) {
    if      (i_face == 0) { ids.resize(3);  ids[0] = pts[2];  ids[1] = pts[1];  ids[2] = pts[0]; }
    else if (i_face == 1) { ids.resize(3);  ids[0] = pts[1];  ids[1] = pts[3];  ids[2] = pts[0]; }
    else if (i_face == 2) { ids.resize(3);  ids[0] = pts[3];  ids[1] = pts[2];  ids[2] = pts[0]; }
    else if (i_face == 3) { ids.resize(3);  ids[0] = pts[2];  ids[1] = pts[3];  ids[2] = pts[1]; }

  } else if (type_cell == VTK_PYRAMID) {
    if      (i_face == 0) { ids.resize(4);  ids[0] = pts[0];  ids[1] = pts[3];  ids[2] = pts[2];  ids[3] = pts[1]; }
    else if (i_face == 1) { ids.resize(3);  ids[0] = pts[0];  ids[1] = pts[1];  ids[2] = pts[4]; }
    else if (i_face == 2) { ids.resize(3);  ids[0] = pts[1];  ids[1] = pts[2];  ids[2] = pts[4]; }
    else if (i_face == 3) { ids.resize(3);  ids[0] = pts[2];  ids[1] = pts[3];  ids[2] = pts[4]; }
    else if (i_face == 4) { ids.resize(3);  ids[0] = pts[3];  ids[1] = pts[0];  ids[2] = pts[4]; }

  } else if (type_cell == VTK_WEDGE) {
    if      (i_face == 0) { ids.resize(3);  ids[0] = pts[0];  ids[1] = pts[1];  ids[2] = pts[2]; }
    else if (i_face == 1) { ids.resize(3);  ids[0] = pts[3];  ids[1] = pts[5];  ids[2] = pts[4]; }
    else if (i_face == 2) { ids.resize(4);  ids[0] = pts[3];  ids[1] = pts[4];  ids[2] = pts[1];  ids[3] = pts[0]; }
    else if (i_face == 3) { ids.resize(4);  ids[0] = pts[1];  ids[1] = pts[4];  ids[2] = pts[5];  ids[3] = pts[2]; }
    else if (i_face == 4) { ids.resize(4);  ids[0] = pts[0];  ids[1] = pts[2];  ids[2] = pts[5];  ids[3] = pts[3]; }

  } else if (type_cell == VTK_HEXAHEDRON) {
    if      (i_face == 0) { ids.resize(4);  ids[0] = pts[0];  ids[1] = pts[3];  ids[2] = pts[2];  ids[3] = pts[1]; }
    else if (i_face == 1) { ids.resize(4);  ids[0] = pts[4];  ids[1] = pts[5];  ids[2] = pts[6];  ids[3] = pts[7]; }
    else if (i_face == 2) { ids.resize(4);  ids[0] = pts[0];  ids[1] = pts[1];  ids[2] = pts[5];  ids[3] = pts[4]; }
    else if (i_face == 3) { ids.resize(4);  ids[0] = pts[3];  ids[1] = pts[7];  ids[2] = pts[6];  ids[3] = pts[2]; }
    else if (i_face == 4) { ids.resize(4);  ids[0] = pts[0];  ids[1] = pts[4];  ids[2] = pts[7];  ids[3] = pts[3]; }
    else if (i_face == 5) { ids.resize(4);  ids[0] = pts[1];  ids[1] = pts[2];  ids[2] = pts[6];  ids[3] = pts[5]; }

  } else if (type_cell == VTK_POLYHEDRON) {
    vtkIdType num_faces = pts[0];
    if (i_face >= num_faces) {
      EG_BUG;
    }
    int id = 1;
    for (int i = 0; i < i_face; ++i) {
      id += pts[id] + 1;
    }
    ids.resize(pts[id]);
    for (int i = 0; i < ids.size(); ++i) {
      ids[i] = pts[id + 1 + i];
    }

  } else {
    EG_BUG; // not implemented
  }
}

vec3_t EgVtkObject::getNormalOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face)
{
  QVector<vtkIdType> ids;
  getFaceOfCell(grid, id_cell, i_face, ids);
  if (ids.size() == 0) {
    EG_BUG;
  }
  QVector<vec3_t> x(ids.size() + 1);
  vec3_t xc(0,0,0);
  for (int i = 0; i < ids.size(); ++i) {
    grid->GetPoint(ids[i], x[i].data());
    xc += x[i];
  }
  x[ids.size()] = x[0];
  xc *= 1.0/ids.size();
  vec3_t n(0,0,0);
  for (int i = 0; i < ids.size(); ++i) {
    vec3_t u = x[i] - xc;
    vec3_t v = x[i+1] - xc;
    n += 0.5*u.cross(v);
  }
  return n;
}

vec3_t EgVtkObject::getCentreOfCellFace(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_face)
{
  QVector<vtkIdType> ids;
  getFaceOfCell(grid, id_cell, i_face, ids);
  if (ids.size() == 0) {
    EG_BUG;
  }
  vec3_t xc(0,0,0), x;
  for (int i = 0; i < ids.size(); ++i) {
    grid->GetPoint(ids[i], x.data());
    xc += x;
  }
  xc *= 1.0/ids.size();
  return xc;
}

void EgVtkObject::getEdgeOfCell(vtkUnstructuredGrid *grid, vtkIdType id_cell, int i_edge, QVector<vtkIdType> &ids)
{
  vtkIdType type_cell = grid->GetCellType(id_cell);
  ids.clear();
  vtkIdType *pts, N_pts;
  grid->GetCellPoints(id_cell, N_pts, pts);
  if (type_cell == VTK_TETRA) {
    ids.resize(2);
    if      (i_edge == 0) { ids[0] = pts[0]; ids[1] = pts[1]; }
    else if (i_edge == 1) { ids[0] = pts[0]; ids[1] = pts[2]; }
    else if (i_edge == 2) { ids[0] = pts[0]; ids[1] = pts[3]; }
    else if (i_edge == 3) { ids[0] = pts[1]; ids[1] = pts[2]; }
    else if (i_edge == 4) { ids[0] = pts[1]; ids[1] = pts[3]; }
    else if (i_edge == 5) { ids[0] = pts[2]; ids[1] = pts[3]; }
  } else {
    EG_BUG; // not implemented
  }
}

bool EgVtkObject::saveGrid(vtkUnstructuredGrid* a_grid, QString file_name)
{
  file_name += ".vtu";
  addVtkTypeInfo(a_grid);
  createIndices(a_grid);
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
  vtu->SetFileName(qPrintable(file_name));
  vtu->SetDataModeToBinary();
  vtu->SetInput(a_grid);
  vtu->Write();
  if(vtu->GetErrorCode()) {
    return false;
  }
  else {
    return true;
  }
}

void EgVtkObject::addVtkTypeInfo(vtkUnstructuredGrid* a_grid)
{
  EG_VTKSP(vtkIntArray, vtk_type);
  vtk_type->SetName("vtk_type");
  vtk_type->SetNumberOfValues(a_grid->GetNumberOfCells());
  //EG_VTKDCC(vtkDoubleArray, cell_VA, a_grid, "cell_VA");
  for (vtkIdType cellId = 0; cellId < a_grid->GetNumberOfCells(); ++cellId) {
    vtk_type->SetValue(cellId, a_grid->GetCellType(cellId));
    //cell_VA->SetValue(cellId, GeometryTools::cellVA(a_grid, cellId, true));
  }
  a_grid->GetCellData()->AddArray(vtk_type);
}

vtkIdType EgVtkObject::addGrid(vtkUnstructuredGrid *main_grid, vtkUnstructuredGrid *grid_to_add, vtkIdType offset)
{
  for (vtkIdType id_node = 0; id_node < grid_to_add->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    grid_to_add->GetPoints()->GetPoint(id_node, x.data());
    main_grid->GetPoints()->SetPoint(offset + id_node, x.data());
    copyNodeData(grid_to_add, id_node, main_grid, offset + id_node);
  }
  for (vtkIdType id_cell = 0; id_cell < grid_to_add->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *pts;
    vtkIdType type_cell = grid_to_add->GetCellType(id_cell);
    grid_to_add->GetCellPoints(id_cell, N_pts, pts);
    QVector <vtkIdType> new_pts(N_pts);
    for(int i=0;i<N_pts;i++) new_pts[i] = offset + pts[i];
    vtkIdType id_new_cell = main_grid->InsertNextCell(type_cell, N_pts, new_pts.data());
    copyCellData(grid_to_add, id_cell, main_grid, id_new_cell);
  }
  return( offset + grid_to_add->GetNumberOfPoints() );
}

QSet<int> EgVtkObject::getAllBoundaryCodes(vtkUnstructuredGrid *grid)
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  QSet<int> bcs;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      bcs.insert(cell_code->GetValue(id_cell));
    }
  }
  return bcs;
}

bool EgVtkObject::cellContainsNode(vtkUnstructuredGrid *grid, vtkIdType id_cell, vtkIdType id_node)
{
  /*
  vtkIdType num_pts, *pts;
  grid->GetCellPoints(id_cell, num_pts, pts);
  for (int i = 0; i < num_pts; ++i) {
    if (pts[i] == id_node) {
      return true;
    }
  }
  return false;
  */

  QVector<vtkIdType> pts;
  getPointsOfCell(grid, id_cell, pts);
  if (pts.contains(id_node)) {
    return true;
  }
  return false;

}

void EgVtkObject::createPolyDataC2C(vtkPolyData *poly_data, QVector<QVector<vtkIdType> > &c2c)
{
  c2c.clear();
  c2c.resize(poly_data->GetNumberOfCells());
  QVector<QSet<vtkIdType> > n2c;
  createPolyDataN2C(poly_data, n2c);
  EG_FORALL_CELLS(id_poly, poly_data) {
    EG_GET_CELL(id_poly, poly_data);
    QSet<vtkIdType> ids;
    for (int i = 0; i < num_pts; ++i) {
      vtkIdType p1 = pts[i];
      vtkIdType p2 = pts[0];
      if (i < num_pts - 1) {
        p2 = pts[i+1];
      }
      QSet<vtkIdType> shared_polys = n2c[p1];
      shared_polys.intersect(n2c[p2]);
      if (shared_polys.size() == 0) {
        EG_BUG;
      }
      if (shared_polys.size() > 2) {
        EG_BUG;
      }
      foreach (vtkIdType id_neigh, shared_polys) {
        if (id_neigh != id_poly) {
          ids.insert(id_neigh);
        }
      }
    }
    {
      c2c[id_poly].resize(ids.size());
      int i = 0;
      foreach (vtkIdType id_neigh, ids) {
        c2c[id_poly][i] = id_neigh;
        ++i;
      }
    }
  }
}

void EgVtkObject::createPolyDataN2C(vtkPolyData *poly_data, QVector<QSet<vtkIdType> > &n2c)
{
  n2c.clear();
  n2c.resize(poly_data->GetNumberOfPoints());
  EG_FORALL_CELLS(id_face, poly_data) {
    EG_GET_CELL(id_face, poly_data);
    for (int i = 0; i < num_pts; ++i) {
      n2c[pts[i]].insert(id_face);
    }
  }
}

void EgVtkObject::createPolyDataN2N(vtkPolyData *poly_data, QVector<QSet<vtkIdType> > &n2n)
{
  n2n.clear();
  n2n.resize(poly_data->GetNumberOfPoints());
  EG_FORALL_CELLS(id_face, poly_data) {
    EG_GET_CELL(id_face, poly_data);
    n2n[pts[0]].insert(pts[num_pts - 1]);
    n2n[pts[num_pts-1]].insert(pts[0]);
    for (int i = 1; i < num_pts; ++i) {
      n2n[pts[i]].insert(pts[i-1]);
      n2n[pts[i-1]].insert(pts[i]);
    }
  }
}

void EgVtkObject::checkGridConsitency(vtkUnstructuredGrid *grid)
{
  vtkIdType num_cells = grid->GetNumberOfCells();
  vtkIdType num_nodes = grid->GetNumberOfPoints();
  for (vtkIdType id_cell = 0; id_cell < num_cells; ++id_cell) {
    QList<vtkIdType> pts;
    getPointsOfCell(grid, id_cell, pts);
    foreach (vtkIdType id_node, pts) {
      if (id_node >= num_nodes) {
        EG_BUG;
      }
    }
  }
}









