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
  m_FeatureAngle = 180.0;
  m_NewBoundaryCode = 1;
  m_OldBoundaryCode = -1;
  setSurfaceIteration();
  setQuickSave(true);
}

void SetBoundaryCode::pass1()
{
  if(!(m_SelectAllVisible || m_OnlyPickedCell || m_OnlyPickedCellAndNeighbours)) {
    using namespace GeometryTools;
    double fa = m_FeatureAngle*M_PI/180.0;
  
    QSet <int> DBC;
    GuiMainWindow::pointer()->getDisplayBoundaryCodes(DBC);
    DBC.insert(9999);
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    
    for (int i = 0; i < pair.size(); ++i) {
      int bc2 = cell_code->GetValue(pair[i].item2);
      if (m_ProcessAll) {
        if (bc2 == m_OldBoundaryCode || m_OldBoundaryCode < 0) {
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
      } else {
        if (DBC.contains(bc2)) {
          if (bc2 == m_OldBoundaryCode || m_OldBoundaryCode < 0) {
            vec3_t n1 = cellNormal(m_Grid, pair[i].item1);
            vec3_t n2 = cellNormal(m_Grid, pair[i].item2);
            n1.normalise();
            n2.normalise();
            if (GeometryTools::angle(n1, n2) > fa) {
              pair[i].terminate = true;
            } else {
              pair[i].terminate = false;
            }
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
  if(m_SelectAllVisible) {
    QSet <int> DBC;
    GuiMainWindow::pointer()->getDisplayBoundaryCodes(DBC);
    DBC.insert(m_NewBoundaryCode);
    foreach(cellId, cells) {
      int bc = cell_code->GetValue(cellId);
      if(DBC.contains(bc)){
        cell_code->SetValue(cellId, m_NewBoundaryCode);
      }
    }
  } else if (m_OnlyPickedCell) {
    cout<<"this->getStart()="<<this->getStart()<<endl;
    cell_code->SetValue(this->getStart(), m_NewBoundaryCode);
  } else if (m_OnlyPickedCellAndNeighbours) {
    cell_code->SetValue(this->getStart(), m_NewBoundaryCode);
    foreach (cellId, c2c[this->getStart()]) {
      cell_code->SetValue(cellId, m_NewBoundaryCode);
    }
  } else {
    foreach (cellId, item) {
      cell_code->SetValue(cellId, m_NewBoundaryCode);
    }
  }
  m_Grid->Modified();
}

