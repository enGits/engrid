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

#include "surfacenodemovementcheck.h"
#include "geometrytools.h"

#include "vtkDelaunay3D.h"
#include "vtkXMLUnstructuredGridWriter.h"

SurfaceNodeMovementCheck::SurfaceNodeMovementCheck(vec3_t x)
{
  m_AuxGrid = vtkUnstructuredGrid::New();
  m_UpdateRequired = true;
  m_OldCentre = x;
}

SurfaceNodeMovementCheck::~SurfaceNodeMovementCheck()
{
  m_AuxGrid->Delete();
}

void SurfaceNodeMovementCheck::addNode(vec3_t x)
{
  m_UpdateRequired = true;
  m_Nodes.append(x);
}

void SurfaceNodeMovementCheck::update()
{
  EG_VTKSP(vtkUnstructuredGrid, grid);
  allocateGrid(grid, 1, m_Nodes.size() + 1, false);
  grid->GetPoints()->SetPoint(0, m_OldCentre.data());
  {
    vtkIdType id_node = 1;
    foreach (vec3_t x, m_Nodes) {
      grid->GetPoints()->SetPoint(id_node, x.data());
      ++id_node;
    }
  }
  EG_VTKSP(vtkDelaunay3D, delaunay);
  delaunay->SetInput(grid);
  delaunay->Update();
  makeCopy(delaunay->GetOutput(), m_AuxGrid, false);
  m_UpdateRequired = false;
}

void SurfaceNodeMovementCheck::write(QString file_name)
{
  EG_VTKSP(vtkXMLUnstructuredGridWriter, vtu);
  vtu->SetFileName(qPrintable(file_name + ".vtu"));
  vtu->SetInput(m_AuxGrid);
  vtu->Write();
}

bool SurfaceNodeMovementCheck::operator()(vec3_t x)
{
  if (m_UpdateRequired) {
    update();
  }
  m_AuxGrid->GetPoints()->SetPoint(0, x.data());
  EG_FORALL_CELLS(id_cell, m_AuxGrid) {
    if (GeometryTools::cellVA(m_AuxGrid, id_cell) < 0) {
      return false;
    }
  }
  return true;
}
