//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "guinormalextrusion.h"
#include "vtkEgNormalExtrusion.h"
#include "containertricks.h"

void GuiNormalExtrusion::before()
{
  populateBoundaryCodes(ui.listWidget);
}

void GuiNormalExtrusion::operate()
{
  EG_VTKSP(vtkEgNormalExtrusion, extr);
  QVector<double> y;
  
  if (ui.radioButtonSimple->isChecked()) {
    y.resize(ui.lineEditSimpleNumLayers->text().toInt() + 1);
    double h = ui.lineEditSimpleHeight->text().toDouble();
    double f = ui.lineEditSimpleIncrease->text().toDouble();
    y[0] = 0.0;
    for (int i = 1; i < y.size(); ++i) {
      y[i] = y[i-1] + h;
      h *= f;
    }
  } else if (ui.radioButtonFixedHeights->isChecked()) {
    y.resize(ui.lineEditFixedHeightsNumLayers->text().toInt() + 1);
    QVector<double> x(y.size());
    for (int i = 0; i < x.size(); ++i) {
      x[i] = i*1.0/(x.size() - 1);
    }
    mat3_t A;
    clinit(A[0]) = pow(x[1],5.0), pow(x[1],3.0), x[1];
    clinit(A[1]) = pow(x[x.size() - 2],5.0), pow(x[x.size() - 2],3.0), x[x.size() - 2];
    clinit(A[2]) = pow(x[x.size() - 1],5.0), pow(x[x.size() - 1],3.0), x[x.size() - 1];
    vec3_t h;
    h[0] = ui.lineEditFixedHeightsHeightFirst->text().toDouble();
    h[2] = ui.lineEditFixedHeightsTotalHeight->text().toDouble();
    h[1] = h[2] - ui.lineEditFixedHeightsHeightLast->text().toDouble();
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
  
  if (ui.radioButtonFixed->isChecked()) {
    extr->SetNormal(vec3_t(ui.lineEditFixedNX->text().toDouble(),
                           ui.lineEditFixedNY->text().toDouble(),
                           ui.lineEditFixedNZ->text().toDouble()));
    double min_dist = ui.lineEditFixedDist->text().toDouble();
    if (min_dist <= 0) {
      extr->SetFixed();
    } else {
      extr->SetPlanar();
      extr->SetMinDist(min_dist);
    }
  }
  if (ui.radioButtonCylinder->isChecked()) {
    extr->SetCylindrical();
    extr->SetOrigin(vec3_t(ui.lineEditCylinderX0->text().toDouble(),
                           ui.lineEditCylinderY0->text().toDouble(),
                           ui.lineEditCylinderZ0->text().toDouble()));
    extr->SetAxis(vec3_t(ui.lineEditCylinderNX->text().toDouble(),
                         ui.lineEditCylinderNY->text().toDouble(),
                         ui.lineEditCylinderNZ->text().toDouble()));
  }
  if (ui.radioButtonRotation->isChecked()) {
    extr->SetRotation();
    extr->SetOrigin(vec3_t(ui.lineEditCylinderX0->text().toDouble(),
                           ui.lineEditCylinderY0->text().toDouble(),
                           ui.lineEditCylinderZ0->text().toDouble()));
    extr->SetAxis(vec3_t(ui.lineEditCylinderNX->text().toDouble(),
                         ui.lineEditCylinderNY->text().toDouble(),
                         ui.lineEditCylinderNZ->text().toDouble()));
  }
  
  QSet<int> bcs;
  getSelectedItems(ui.listWidget, bcs);
  extr->SetBoundaryCodes(bcs);
  EG_VTKSP(vtkUnstructuredGrid,ug);
  makeCopy(grid, ug);
  extr->SetInput(ug);
  extr->Update();
  makeCopy(extr->GetOutput(), grid);
}
