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
#include "setboundarycode.h"
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include "guimainwindow.h"

SetBoundaryCode::SetBoundaryCode()
{
  feature_angle = 180.0;
  boundary_code = 1;
  setSurfaceIteration();
  setQuickSave(true);
}

void SetBoundaryCode::pass1()
{
  if(!(SelectAllVisible || OnlyPickedCell || OnlyPickedCellAndNeighbours)) {
    using namespace GeometryTools;
    double fa = feature_angle*M_PI/180.0;
  
    QSet <int> DBC;
    GuiMainWindow::pointer()->getDisplayBoundaryCodes(DBC);
    DBC.insert(9999);
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    
    for (int i = 0; i < pair.size(); ++i) {
      int bc1 = cell_code->GetValue(pair[i].item1);
      int bc2 = cell_code->GetValue(pair[i].item2);
      if (ProcessAll) {
        vec3_t n1 = cellNormal(m_Grid, pair[i].item1);
        vec3_t n2 = cellNormal(m_Grid, pair[i].item2);
        n1.normalise();
        n2.normalise();
        if (GeometryTools::angle(n1, n2) > fa) {
          pair[i].terminate = true;
        } else {
          pair[i].terminate = false;
        }
      } else {
        if (DBC.contains(bc2)) {
          vec3_t n1 = cellNormal(m_Grid, pair[i].item1);
          vec3_t n2 = cellNormal(m_Grid, pair[i].item2);
          n1.normalise();
          n2.normalise();
          if (GeometryTools::angle(n1, n2) > fa) {
            pair[i].terminate = true;
          } else {
            pair[i].terminate = false;
          }
        } else {
          pair[i].terminate = true;
        }
      }
    }
  }
}

void SetBoundaryCode::pass2()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  vtkIdType cellId;
  if(SelectAllVisible) {
    QSet <int> DBC;
    GuiMainWindow::pointer()->getDisplayBoundaryCodes(DBC);
    DBC.insert(boundary_code);
    foreach(cellId, cells) {
      int bc = cell_code->GetValue(cellId);
      if(DBC.contains(bc)){
        cell_code->SetValue(cellId, boundary_code);
      }
    }
  } else if (OnlyPickedCell) {
    cout<<"this->getStart()="<<this->getStart()<<endl;
    cell_code->SetValue(this->getStart(), boundary_code);
  } else if (OnlyPickedCellAndNeighbours) {
    cell_code->SetValue(this->getStart(), boundary_code);
    foreach (cellId, c2c[this->getStart()]) {
      cell_code->SetValue(cellId, boundary_code);
    }
  } else {
    foreach (cellId, item) {
      cell_code->SetValue(cellId, boundary_code);
    }
  }
  m_Grid->Modified();
}

