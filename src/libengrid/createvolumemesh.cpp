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

#include "createvolumemesh.h"
#include "deletetetras.h"
#include "guimainwindow.h"
#include "updatedesiredmeshdensity.h"
#include "createboundarylayer.h"

#include <vtkXMLUnstructuredGridWriter.h>

CreateVolumeMesh::CreateVolumeMesh()
{
  EG_TYPENAME;
}

int CreateVolumeMesh::numVolumeCells()
{
  int N = 0;
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isVolume(id_cell, m_Grid)) {
      ++N;
    }
  }
  return N;
}

void CreateVolumeMesh::operate()
{
  readSettings();
  m_Part.trackGrid(m_Grid);
  double a = m_MaximalEdgeLength;
  double V = a*a*a/(6*sqrt(2.0));
  int N1 = 0;
  int N2 = numVolumeCells();
  int pass = 1;
  QString V_txt;
  V_txt.setNum(V);
  bool done = false;
  while (!done) {
    N1 = N2;
    QString flags;
    QString q_txt = "10/0";
    if (pass > 1) {
      q_txt = "1.2/0";
    }
    flags = QString("pq") + q_txt + "a" + V_txt;
    if (N2 > 0) {
      flags += "m";
    }
    cout << "TetGen pass " << pass << "  flags=" << qPrintable(flags) << endl;
    tetgen(flags);
    N2 = numVolumeCells();
    cout << N2 << endl;
    //if (pass > 1) N2 = N1;
    ++pass;
    if (fabs(double(N2-N1)/N1) < 0.05 || pass > 1) {
      done = true;
    }
  }
  CreateBoundaryLayer blayer;
  blayer.setGrid(m_Grid);
  blayer.setAllCells();
  blayer();
}

