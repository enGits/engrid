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
#include "guinormalextrusion.h"
#include "vtkEgNormalExtrusion.h"
#include "containertricks.h"
#include "guimainwindow.h"

void GuiNormalExtrusion::before()
{
  populateBoundaryCodes(m_Ui.listWidget);
}

void GuiNormalExtrusion::operate()
{
  QSet<int> bcs;
  getSelectedItems(m_Ui.listWidget, bcs);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  QSet<int> volume_codes;
  EG_FORALL_CELLS(id_cell, m_Grid) {
    if (isSurface(id_cell, m_Grid)) {
      if (bcs.contains(cell_code->GetValue(id_cell))) {
        vtkIdType id_volume_cell = m_Part.getVolumeCell(id_cell);
        if (id_volume_cell != -1) {
          volume_codes.insert(cell_code->GetValue(id_volume_cell));
        }
      }
    }
  }

  EG_VTKSP(vtkEgNormalExtrusion, extr);
  QVector<double> y;
  
  if (m_Ui.radioButtonSimple->isChecked()) {
    y.resize(m_Ui.lineEditSimpleNumLayers->text().toInt() + 1);
    double h = m_Ui.lineEditSimpleHeight->text().toDouble();
    double f = m_Ui.lineEditSimpleIncrease->text().toDouble();
    y[0] = 0.0;
    for (int i = 1; i < y.size(); ++i) {
      y[i] = y[i-1] + h;
      h *= f;
    }
  } else if (m_Ui.radioButtonFixedHeights->isChecked()) {
    y.resize(m_Ui.lineEditFixedHeightsNumLayers->text().toInt() + 1);
    QVector<double> x(y.size());
    for (int i = 0; i < x.size(); ++i) {
      x[i] = i*1.0/(x.size() - 1);
    }
    mat3_t A;
    clinit(A[0]) = pow(x[1],5.0), pow(x[1],3.0), x[1];
    clinit(A[1]) = pow(x[x.size() - 2],5.0), pow(x[x.size() - 2],3.0), x[x.size() - 2];
    clinit(A[2]) = pow(x[x.size() - 1],5.0), pow(x[x.size() - 1],3.0), x[x.size() - 1];
    vec3_t h;
    h[0] = m_Ui.lineEditFixedHeightsHeightFirst->text().toDouble();
    h[2] = m_Ui.lineEditFixedHeightsTotalHeight->text().toDouble();
    h[1] = h[2] - m_Ui.lineEditFixedHeightsHeightLast->text().toDouble();
    mat3_t AI = A.inverse();
    vec3_t coeff = AI*h;
    for (int i = 0; i < y.size(); ++i) {
      y[i] = coeff[0]*pow(x[i],5.0) + coeff[1]*pow(x[i],3.0) + coeff[2]*x[i];
      if (i > 0) {
        if (y[i] < y[i-1]) {
          EG_ERR_RETURN("unable to compute layer heights");
        }
      }
    }
  }
  extr->SetLayers(y);
  
  if (m_Ui.radioButtonFixed->isChecked()) {
    extr->SetNormal(vec3_t(m_Ui.lineEditFixedNX->text().toDouble(),
                           m_Ui.lineEditFixedNY->text().toDouble(),
                           m_Ui.lineEditFixedNZ->text().toDouble()));
    double min_dist = m_Ui.lineEditFixedDist->text().toDouble();
    if (min_dist <= 0) {
      extr->SetFixed();
    } else {
      extr->SetPlanar();
      extr->SetMinDist(min_dist);
    }
  }
  if (m_Ui.radioButtonCylinder->isChecked()) {
    extr->SetCylindrical();
    extr->SetOrigin(vec3_t(m_Ui.lineEditCylinderX0->text().toDouble(),
                           m_Ui.lineEditCylinderY0->text().toDouble(),
                           m_Ui.lineEditCylinderZ0->text().toDouble()));
    extr->SetAxis(vec3_t(m_Ui.lineEditCylinderNX->text().toDouble(),
                         m_Ui.lineEditCylinderNY->text().toDouble(),
                         m_Ui.lineEditCylinderNZ->text().toDouble()));
  }
  if (m_Ui.radioButtonRotation->isChecked()) {
    extr->SetRotation();
    extr->SetOrigin(vec3_t(m_Ui.lineEditCylinderX0->text().toDouble(),
                           m_Ui.lineEditCylinderY0->text().toDouble(),
                           m_Ui.lineEditCylinderZ0->text().toDouble()));
    extr->SetAxis(vec3_t(m_Ui.lineEditCylinderNX->text().toDouble(),
                         m_Ui.lineEditCylinderNY->text().toDouble(),
                         m_Ui.lineEditCylinderNZ->text().toDouble()));
  }

  if (m_Ui.radioButtonNoRestrict->isChecked()) {
    extr->SetRestrictNone();
  }
  if (m_Ui.radioButtonXY->isChecked()) {
    extr->SetRestrictXY();
  }
  if (m_Ui.radioButtonXZ->isChecked()) {
    extr->SetRestrictXZ();
  }
  if (m_Ui.radioButtonYZ->isChecked()) {
    extr->SetRestrictYZ();
  }

  if (m_Ui.checkBoxNewVolume->isChecked()) {
    extr->SetRemoveInternalFacesOff();
  }

  QList<VolumeDefinition> vols = GuiMainWindow::pointer()->getAllVols();
  int vc = 1;
  foreach (VolumeDefinition vol, vols) {
    vc = max(vc, vol.getVC());
  }
  if (m_Ui.checkBoxNewVolume->isChecked()) {
    EG_FORALL_CELLS(id_cell, m_Grid) {
      if (isVolume(id_cell, m_Grid)) {
        if (cell_code->GetValue(id_cell) == 0) {
          cell_code->SetValue(id_cell, vc);
        }
      }
    }
  }


  QSet<int> old_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  extr->SetBoundaryCodes(bcs);
  EG_VTKSP(vtkUnstructuredGrid,ug);
  makeCopy(m_Grid, ug);
  extr->SetInputData(ug);
  extr->Update();
  makeCopy(extr->GetOutput(), m_Grid);
  QSet<int> new_bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();

  if (m_Ui.checkBoxNewVolume->isChecked()) {
    ++vc;
    VolumeDefinition new_vol(m_Ui.lineEditVolumeName->text(), vc);
    foreach (int bc, new_bcs) {
      if (bcs.contains(bc)) {
        new_vol.addBC(bc, -1);
      } else if (!old_bcs.contains(bc)) {
        new_vol.addBC(bc, 1);
      }
    }
    vols.append(new_vol);
    GuiMainWindow::pointer()->setAllVols(vols);
    EG_FORALL_CELLS(id_cell, m_Grid) {
      if (isVolume(id_cell, m_Grid)) {
        if (cell_code->GetValue(id_cell) == 0) {
          cell_code->SetValue(id_cell, vc);
        }
      }
    }
  } else {
    if (volume_codes.size() == 1) {
      EG_FORALL_CELLS(id_cell, m_Grid) {
        if (isVolume(id_cell, m_Grid)) {
          if (cell_code->GetValue(id_cell) == 0) {
            cell_code->SetValue(id_cell, *(volume_codes.begin()));
          }
        }
      }
    }
  }
}
