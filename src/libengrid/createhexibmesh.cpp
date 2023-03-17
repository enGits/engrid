// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "createhexibmesh.h"
#include "engrid.h"
#include "guimainwindow.h"
#include "updatedesiredmeshdensity.h"
#include "deletevolumegrid.h"

CreateHexIbMesh::CreateHexIbMesh()
{
  m_MinDim = 30;
  m_InsidePosition = vec3_t(0, 0, 0);
}

double CreateHexIbMesh::meshSize(vtkIdType id_face)
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired");
  EG_GET_CELL(id_face, m_Grid);
  double h = 0;
  for (int i = 0; i < 3; ++i) {
    h = max(h,cl->GetValue(pts[i]));
  }
  return h;
}

QString CreateHexIbMesh::bigIntText(long long int N)
{
  QString text, num;
  text.setNum(N%1000);
  N /= 1000;
  int L = 3;
  while (N > 0) {
    num.setNum(N%1000);
    text = num + "," + text.rightJustified(L, '0');
    L += 4;
    N /= 1000;
  }
  return text;
}

double CreateHexIbMesh::meshSize(const QList<vtkIdType> &faces)
{
  double mesh_size = m_MaxEdgeLength;
  foreach (vtkIdType id_face, faces) {
    mesh_size = min(meshSize(id_face), mesh_size);
  }
  return mesh_size;
}

int CreateHexIbMesh::refine()
{
  int N = 0;
  m_Octree.resetRefineMarks();
  int old_num_cells = m_Octree.getNumCells();
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (!m_Octree.hasChildren(cell)) {
      double H = m_Octree.getDx(cell);
      if (H > 2*m_MeshSize[cell]*m_MinDim && H > 2*m_MinSize) {
        m_Octree.markToRefine(cell);
        ++N;
      }
    }
  }
  m_Octree.refineAll();
  m_Faces.insert(old_num_cells, m_Octree.getNumCells() - old_num_cells, QList<vtkIdType>());
  m_MeshSize.insert(old_num_cells, m_Octree.getNumCells() - old_num_cells, m_MaxEdgeLength);
  for (int cell = old_num_cells; cell < m_Octree.getNumCells(); ++cell) {
    foreach (vtkIdType id_face, m_Faces[m_Octree.getParent(cell)]) {
      EG_GET_CELL(id_face, m_Grid);
      if (num_pts != 3) {
        EG_BUG;
      }
      QVector<vec3_t> tri(3);
      for (int i = 0; i < 3; ++i) {
        m_Grid->GetPoint(pts[i], tri[i].data());
      }

      // compute scale factor to make sure that we have the required minimal number of layers
      double H = max(m_Octree.getDx(cell), max(m_Octree.getDy(cell), m_Octree.getDz(cell)));
      double h = meshSize(id_face)*m_MinNumLayersWithRequiredResolution;
      double scale = 1.0;
      if (h < H) {
        scale += h/H;
      }
      if (m_Octree.triangleIntersectsCell(cell, tri, scale)) {
        m_Faces[cell].append(id_face);
      }
    }
    m_MeshSize[cell] = min(meshSize(m_Faces[cell]), m_ELSManager.minEdgeLength(m_Octree.getCellCentre(cell)));
  }
  int num_supercells = 0;
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (!m_Octree.hasChildren(cell)) {
      ++num_supercells;
    }
  }
  cout << qPrintable(bigIntText(num_supercells)) << " super cells" << endl;
  return N;
}

void CreateHexIbMesh::updateMeshSize()
{
  // find smallest cell size (on highest level)
  double H = 1e99;
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (!m_Octree.hasChildren(cell)) {
      double l = m_Octree.getLevel(cell);
      H = min(H, m_MeshSize[cell]*pow(2.0, l));
    }
  }

  // set discrete cell size according to highest level
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (!m_Octree.hasChildren(cell)) {
      double l = m_Octree.getLevel(cell);
      m_MeshSize[cell] = H/pow(2.0, l);
    }
  }
}

void CreateHexIbMesh::findInsideCells(MeshPartition &part, QList<vtkIdType> &inside_cells)
{
  vtkUnstructuredGrid *grid = part.getGrid();
  QVector<bool> is_inside(grid->GetNumberOfCells(), false);

  // recompute faces for cells with no children with a minimal scale factor
  for (int cell = 0; cell < m_Octree.getNumCells(); ++cell) {
    if (m_Octree.getLevel(cell) > 0 && !m_Octree.hasChildren(cell)) {
      m_Faces[cell].clear();
      foreach (vtkIdType id_face, m_Faces[m_Octree.getParent(cell)]) {
        EG_GET_CELL(id_face, m_Grid);
        QVector<vec3_t> tri(3);
        for (int i = 0; i < 3; ++i) {
          m_Grid->GetPoint(pts[i], tri[i].data());
        }
        if (m_Octree.triangleIntersectsCell(cell, tri, 1.01)) {
          m_Faces[cell].append(id_face);
        }
      }
    }
  }

  // find closest cell to m_InsidePosition and mark boundary cells
  vtkIdType id_start = -1;
  {
    double min_dist = 1e99;
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      vec3_t x = cellCentre(grid, id_cell);
      double dist = (m_InsidePosition - x).abs();
      if (dist < min_dist) {
        min_dist = dist;
        id_start = id_cell;
      }
      int cell = m_Octree.findCell(x);
      if (m_Faces[cell].size() > 0) {
        is_inside[id_cell] = true;
      }
    }
  }
  if (id_start == -1) {
    EG_BUG;
  }
  QList<vtkIdType> front_cells;
  front_cells << id_start;
  is_inside[id_start] = true;
  while (front_cells.size() != 0) {
    QList<vtkIdType> new_front_cells;
    foreach (vtkIdType id_cell, front_cells) {
      for (int i = 0; i < part.c2cGSize(id_cell); ++i) {
        vtkIdType id_neigh = part.c2cGG(id_cell, i);
        if (id_neigh >= 0) {
          if (!is_inside[id_neigh]) {
            is_inside[id_neigh] = true;
            new_front_cells << id_neigh;
          }
        }
      }
    }
    front_cells = new_front_cells;
  }
  inside_cells.clear();
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (is_inside[id_cell]) {
      inside_cells << id_cell;
    }
  }
}

void CreateHexIbMesh::operate()
{
  DeleteVolumeGrid del_vol;
  del_vol();
  m_Part.setAllCells();
  cout << "refining grid" << endl;
  updateNodeInfo();
  UpdateDesiredMeshDensity update_mesh_density;
  update_mesh_density.readSettings();
  update_mesh_density.setRelaxationOff();
  update_mesh_density();

  QString buffer = GuiMainWindow::pointer()->getXmlSection("engrid/surface/settings").replace("\n", " ");
  if (!buffer.isEmpty()) {
    QTextStream in(&buffer, QIODevice::ReadOnly);
    in >> m_MaxEdgeLength;
    in >> m_MinEdgeLength;
    in >> m_GrowthFactor;
  } else {
    m_MaxEdgeLength = 1000.0;
    m_MinEdgeLength = 0.0;
    m_GrowthFactor = 1.5;
  }
  m_ELSManager.read();

  m_Octree.setSmoothTransitionOn();
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      EG_BUG;
    }
  }
  {
    double bounds[6];
    m_Grid->GetBounds(bounds);
    vec3_t x1(bounds[0], bounds[2], bounds[4]);
    vec3_t x2(bounds[1], bounds[3], bounds[5]);
    vec3_t xc = 0.5*(x1 + x2);
    vec3_t Dx = x2 - xc;
    double dx_max = max(Dx[0], max(Dx[1], Dx[2]));
    Dx = vec3_t(dx_max, dx_max, dx_max);
    x1 = xc - 2*Dx;
    x2 = xc + 2*Dx;
    m_Octree.setBounds(x1, x2);
    m_MinSize = 0;
  }
  m_Faces.resize(1);
  m_MeshSize.resize(1);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    m_Faces[0].append(id_cell);
  }
  m_MeshSize[0] = meshSize(m_Faces[0]);

  int N;
  do {
    N = refine();
  } while (N > 0);
  updateMeshSize();

  cout << "deleting outside super cells" << endl;
  EG_VTKSP(vtkUnstructuredGrid, otgrid);
  m_Octree.toVtkGridPolyhedral(otgrid, true);
  MeshPartition otpart(otgrid);
  otpart.setAllCells();
  QList<vtkIdType> add_cells;
  findInsideCells(otpart, add_cells);
  otpart.setCells(add_cells);
  EG_VTKDCC(vtkDoubleArray, cl, otgrid, "cell_subres");
  double fvcells = 0;
  int N_min = 1000000;
  int N_max = 0;
  foreach (vtkIdType id_cell, add_cells) {
    vec3_t x = cellCentre(otgrid, id_cell);
    int cell = m_Octree.findCell(x);
    double h = m_MeshSize[cell];
    cl->SetValue(id_cell, h);
    double Ni = m_Octree.getDx(cell)/h;
    double Nj = m_Octree.getDy(cell)/h;
    double Nk = m_Octree.getDx(cell)/h;
    fvcells += Ni*Nj*Nk;
    N_min = min(N_min, min(int(Ni), min(int(Nj), int(Nk))));
    N_max = max(N_max, max(int(Ni), max(int(Nj), int(Nk))));
  }
  cout << add_cells.size() << " super cells" << endl;
  if (fvcells > 1e12) {
    cout << fvcells << " finite volume cells!!!!" << endl;
  } else {
    long long int fvcells_li = fvcells;
    cout << qPrintable(bigIntText(fvcells_li)) << " finite volume cells" << endl;
  }
  cout << "smallest dimension: " << N_min << endl;
  cout << "largest dimension: " << N_max << endl;
  m_Part.addPartition(otpart);
  m_Part.setAllCells();
}
