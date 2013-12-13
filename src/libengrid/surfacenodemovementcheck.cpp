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

#include "surfacenodemovementcheck.h"
#include "geometrytools.h"

#include "vtkDelaunay3D.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkKdTreePointLocator.h"

SurfaceNodeMovementCheck::SurfaceNodeMovementCheck()
{
  m_AuxGrid = vtkUnstructuredGrid::New();
  m_Grid = NULL;
}

SurfaceNodeMovementCheck::~SurfaceNodeMovementCheck()
{
  m_AuxGrid->Delete();
}

void SurfaceNodeMovementCheck::setGrid(vtkUnstructuredGrid *grid)
{
  m_Grid = grid;
  EG_VTKSP(vtkUnstructuredGrid, tmp_grid);
  allocateGrid(tmp_grid, 1, m_Grid->GetNumberOfPoints(), false);
  {
    vec3_t x;
    EG_FORALL_NODES(id_node, m_Grid) {
      m_Grid->GetPoint(id_node, x.data());
      tmp_grid->GetPoints()->SetPoint(id_node, x.data());
    }
  }
  EG_VTKSP(vtkDelaunay3D, delaunay);
  delaunay->SetInputData(tmp_grid);
  delaunay->Update();
  makeCopy(delaunay->GetOutput(), m_AuxGrid, false);
  int N1 = m_Grid->GetNumberOfPoints();
  int N2 = m_AuxGrid->GetNumberOfPoints();
  if (N1 != N2) {
    EG_BUG;
  }
  EG_VTKSP(vtkKdTreePointLocator, loc);
  loc->SetDataSet(m_AuxGrid);
  loc->BuildLocator();
  m_IdMap.resize(m_Grid->GetNumberOfPoints());
  EG_FORALL_NODES(id_node, m_Grid) {
    vec3_t x;
    m_Grid->GetPoint(id_node, x.data());
    vtkIdType id_aux_node = id_node;//loc->FindClosestPoint(x.data());
    if (id_aux_node < 0) {
      EG_BUG;
    }
    m_IdMap[id_node] = id_aux_node;
  }
  m_AuxPart.setGrid(m_AuxGrid);
  m_AuxPart.setAllCells();
  /*
  QVector<vtkIdType> nodes(m_AuxGrid->GetNumberOfPoints());
  EG_FORALL_NODES(id_aux_node, m_AuxGrid) {
    nodes[id_aux_node] = id_aux_node;
  }
  m_AuxPart.setNodes(nodes);
  */
  int N3 = 0;
  EG_FORALL_NODES(id_aux_node, m_AuxGrid) {
    if (m_AuxPart.localNode(id_aux_node) == -1) {
      ++N3;
    }
  }
  cout << N3 << endl;
}

bool SurfaceNodeMovementCheck::checkNode(vtkIdType id_node, vec3_t x)
{
  bool good = true;
  vec3_t x_old = x;
  vtkIdType id_aux_node = m_IdMap[id_node];
  if (m_AuxPart.localNode(id_aux_node) == -1) {
    good = false;
  } else {
    m_AuxGrid->GetPoints()->SetPoint(id_aux_node, x.data());
    for (int i = 0; i < m_AuxPart.n2cGSize(id_aux_node); ++i) {
      vtkIdType id_aux_cell = m_AuxPart.n2cGG(id_aux_node, i);
      double V = GeometryTools::cellVA(m_AuxGrid, id_aux_cell);
      if (V < 0) {
        good = false;
        break;
      }
    }
  }
  m_AuxGrid->GetPoints()->SetPoint(id_aux_node, x_old.data());
  return good;
}

void SurfaceNodeMovementCheck::write(QString file_name)
{
  EG_VTKSP(vtkXMLUnstructuredGridWriter, vtu);
  vtu->SetFileName(qPrintable(file_name + ".vtu"));
  vtu->SetInputData(m_AuxGrid);
  vtu->Write();
}

