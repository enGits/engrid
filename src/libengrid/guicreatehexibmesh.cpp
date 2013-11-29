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
#include "guicreatehexibmesh.h"
#include "geometrytools.h"

GuiCreateHexIbMesh::GuiCreateHexIbMesh()
{
}

void GuiCreateHexIbMesh::before()
{
  vec3_t x_centre(0,0,0);
  double A = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      double A_cell = GeometryTools::cellVA(m_Grid, id_cell);
      x_centre += A_cell*cellCentre(m_Grid, id_cell);
      A += A_cell;
    }
  }
  x_centre *= 1.0/A;
  setVector(x_centre, m_Ui.m_LineEditCentre);
}

void GuiCreateHexIbMesh::operate()
{
  m_CreateMesh.setMinNumLayersWithRequiredResolution(m_Ui.m_SpinBoxMinNumLayers->value());
  m_CreateMesh.setMinDim(m_Ui.m_SpinBoxMinDim->value());
  vec3_t x_centre = getVector(m_Ui.m_LineEditCentre);
  m_CreateMesh.setInsidePosition(x_centre);
  m_CreateMesh();
}
