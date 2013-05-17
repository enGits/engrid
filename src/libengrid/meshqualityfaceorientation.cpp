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
#include "meshqualityfaceorientation.h"
#include "geometrytools.h"
#include "guimainwindow.h"

void MeshQualityFaceOrientation::operate()
{
  using namespace GeometryTools;
  EG_VTKDCC(vtkDoubleArray, cell_mesh_quality, m_Grid, "cell_mesh_quality");
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid)) {
      vec3_t x_face = cellCentre(m_Grid, id_cell);
      vec3_t n_face = cellNormal(m_Grid, id_cell);
      n_face.normalise();
      CadInterface* cad_interface = GuiMainWindow::pointer()->getCadInterface(cell_code->GetValue(id_cell), true);
      if (cad_interface) {
        //proj->snapNode(x_face, -1);
        cad_interface->project(x_face, n_face);
        if (cad_interface->failed()) {
          cell_mesh_quality->SetValue(id_cell, 1.0);
        } else {
          vec3_t n_surf = cad_interface->getLastNormal();
          double mq = 0.5*(n_surf*n_face + 1);
          cell_mesh_quality->SetValue(id_cell, mq);
        }
      } else {
        cell_mesh_quality->SetValue(id_cell, 1.0);
      }
    } else {
      cell_mesh_quality->SetValue(id_cell, 1.0);
    }
  }
  computeNodesFromCells();
}
