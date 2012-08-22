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

#include "vtkDelaunay3D.h"
#include "geometrytools.h"

SurfaceNodeMovementCheck::SurfaceNodeMovementCheck()
{
  m_VolumeGrid = NULL;
}

void SurfaceNodeMovementCheck::setSurfaceGrid(vtkUnstructuredGrid *surface_grid)
{
  m_SurfaceGrid = surface_grid;
  m_VolumeGrid = vtkUnstructuredGrid::New();
  m_SurfacePart.setGrid(m_SurfaceGrid);
  update();
}

void SurfaceNodeMovementCheck::update()
{
  EG_VTKSP(vtkDelaunay3D, triangulator);
  triangulator->SetInput(m_SurfaceGrid);
  triangulator->Update();
  makeCopy(triangulator->GetOutput(), m_VolumeGrid);
  m_VolumePart.setGrid(m_VolumeGrid);
  m_VolumePart.setAllCells();
  m_SurfacePart.setGrid(m_SurfaceGrid);
  m_SurfacePart.setAllCells();
  int N1 = m_VolumeGrid->GetNumberOfPoints();
  int N2 = m_SurfaceGrid->GetNumberOfPoints();
  if (m_VolumeGrid->GetNumberOfPoints() < m_SurfaceGrid->GetNumberOfPoints()) {
    EG_BUG;
  }
}

SurfaceNodeMovementCheck::~SurfaceNodeMovementCheck()
{
  m_VolumeGrid->Delete();
}

bool SurfaceNodeMovementCheck::moveNode(vtkIdType id_node, vec3_t x_new)
{
  // reject any node which could not be triangulated (orphan)
  if (m_VolumePart.localNode(id_node) < 0) {
    return false;
  }

  // reject any nodes which are neighbours of orphans
  for (int i = 0; i < m_SurfacePart.n2nGSize(id_node); ++i) {
    if (m_VolumePart.localNode(m_SurfacePart.n2nGG(id_node,i)) < 0) {
      return false;
    }
  }

  // check all adjacent volumes for new node position
  vec3_t x_old;
  m_VolumeGrid->GetPoint(id_node, x_old.data());
  m_VolumeGrid->GetPoints()->SetPoint(id_node, x_new.data());
  for (int i = 0; i < m_VolumePart.n2cGSize(id_node); ++i) {
    vtkIdType id_cell = m_VolumePart.n2cGG(id_node, i);
    double V = GeometryTools::cellVA(m_VolumeGrid, id_cell, true);
    if (V <= 0) {
      m_VolumeGrid->GetPoints()->SetPoint(id_node, x_old.data());
      return false;
    }
  }
  m_SurfaceGrid->GetPoints()->SetPoint(id_node, x_new.data());
  return true;
}
