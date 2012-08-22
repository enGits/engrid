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
  cout << "num nodes: " << m_VolumeGrid->GetNumberOfPoints() << endl;
  cout << "num cells: " << m_VolumeGrid->GetNumberOfCells() << endl;
}

SurfaceNodeMovementCheck::~SurfaceNodeMovementCheck()
{
  m_VolumeGrid->Delete();
}
