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
#include "createcadtesselation.h"
#include "guimainwindow.h"

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>
#include <vtkContourFilter.h>
#include <vtkDecimatePro.h>

CreateCadTesselation::CreateCadTesselation()
{
  m_ScanMemory = 0.5;
}

void CreateCadTesselation::scan()
{
  m_Dx = (m_X2[0] - m_X1[0])/(m_Ni-1);
  m_Dy = (m_X2[1] - m_X1[1])/(m_Nj-1);
  m_Dz = (m_X2[2] - m_X1[2])/(m_Nk-1);
  EG_VTKSP(vtkImageData, gdata);
  gdata->SetDimensions(m_Ni, m_Nj, m_Nk);
  EG_VTKSP(vtkFloatArray, g);
  g->SetName("g");
  g->SetNumberOfValues(m_Ni*m_Nj*m_Nk);
  gdata->GetPointData()->AddArray(g);
  for (int i = 0; i < m_Ni; ++i) {
    for (int j = 0; j < m_Nj; ++j) {
      for (int k = 0; k < m_Nk; ++k) {
        g->SetValue(getIdx(i,j,k), 0);
      }
    }
  }
  m_XScan1 = m_X2;
  m_XScan2 = m_X1;

  vec3_t x_in, x_out, n_in, n_out;

  int progress = 0;
  m_GeometryFound = false;

  // scan ij plane -> k direction
  //      xy plane -> z direction
  for (int i = 0; i < m_Ni; ++i) {
    for (int j = 0; j < m_Nj; ++j) {
      vec3_t x = getX(i,j,0);
      while (shootRay(x, vec3_t(0,0,1), x_in, x_out, n_in, n_out)) {
        m_GeometryFound = true;
        m_XScan1[2] = min(m_XScan1[2], x_in[2]);
        m_XScan2[2] = max(m_XScan2[2], x_out[2]);
        int k1 = int((x_in[2] - m_X1[2])/m_Dz) + 1;
        int k2 = int((x_out[2] - m_X1[2])/m_Dz);
        for (int k = max(0,k1); k <= min(m_Nk-1,k2); ++k) {
          g->SetValue(getIdx(i,j,k), 1);
        }
        x = x_out;
        x[2] += 1e-10*m_Dz;
      }
    }
    progress += m_Nj;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  // scan ik plane -> j direction
  //      xz plane -> y direction
  for (int i = 0; i < m_Ni; ++i) {
    for (int k = 0; k < m_Nk; ++k) {
      vec3_t x = getX(i,0,k);
      while (shootRay(x, vec3_t(0,1,0), x_in, x_out, n_in, n_out)) {
        m_GeometryFound = true;
        m_XScan1[1] = min(m_XScan1[1], x_in[1]);
        m_XScan2[1] = max(m_XScan2[1], x_out[1]);
        int j1 = int((x_in[1]-m_X1[1])/m_Dy) + 1;
        int j2 = int((x_out[1]-m_X1[1])/m_Dy);
        for (int j = max(0,j1); j <= min(m_Nj-1,j2); ++j) {
          g->SetValue(getIdx(i,j,k), 1);
        }
        x = x_out;
        x[1] += 1e-10*m_Dy;
      }
    }
    progress += m_Nk;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  // scan jk plane -> i direction
  //      yz plane -> x direction
  for (int j = 0; j < m_Nj; ++j) {
    for (int k = 0; k < m_Nk; ++k) {
      vec3_t x = getX(0,j,k);
      while (shootRay(x, vec3_t(1,0,0), x_in, x_out, n_in, n_out)) {
        m_GeometryFound = true;
        m_XScan1[0] = min(m_XScan1[0], x_in[0]);
        m_XScan2[0] = max(m_XScan2[0], x_out[0]);
        int i1 = int((x_in[0]-m_X1[0])/m_Dx) + 1;
        int i2 = int((x_out[0]-m_X1[0])/m_Dx);
        for (int i = max(0,i1); i <= min(m_Ni-1,i2); ++i) {
          g->SetValue(getIdx(i,j,k), 1);
        }
        x = x_out;
        x[0] += 1e-10*m_Dx;
      }
    }
    progress += m_Nk;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  gdata->GetPointData()->SetActiveScalars("g");
  EG_VTKSP(vtkContourFilter, contour);
  contour->SetInput(gdata);
  contour->SetNumberOfContours(1);
  contour->SetValue(0, 0.5);
  EG_VTKSP(vtkTriangleFilter, tri);
  tri->SetInput(contour->GetOutput());
  EG_VTKSP(vtkDecimatePro, decimate);
  decimate->PreserveTopologyOn();
  decimate->SetTargetReduction(1.0);
  decimate->SetFeatureAngle(GeometryTools::deg2rad(45));
  decimate->SetInput(tri->GetOutput());
  decimate->Update();

  allocateGrid(m_Grid, decimate->GetOutput()->GetNumberOfPolys(), decimate->GetOutput()->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    decimate->GetOutput()->GetPoint(id_node, x.data());
    m_Grid->GetPoints()->SetPoint(id_node, x.data());
  }
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < decimate->GetOutput()->GetNumberOfPolys(); ++id_cell) {
    EG_GET_CELL(id_cell, decimate->GetOutput());
    vtkIdType id_new_cell = m_Grid->InsertNextCell(type_cell, num_pts, pts);
    cell_code->SetValue(id_new_cell, 1);
  }  
}

void CreateCadTesselation::operate()
{
  int N = int(pow(m_ScanMemory/sizeof(float), 1.0/3.0));
  m_Ni = N;
  m_Nj = N;
  m_Nk = N;
  m_X1 = vec3_t(-1e6,-1e6,-1e6);
  m_X2 = vec3_t( 1e6, 1e6, 1e6);
  double new_range, old_range;
  int count = 0;
  do {
    old_range = (m_X1-m_X2).abs();
    QString text;
    text.setNum((m_X1-m_X2).abs());
    text = QString("scanning (range=") + text + ")";
    GuiMainWindow::pointer()->resetProgress(text, 3*N*N);
    scan();
    if (m_GeometryFound) {
      m_X1 = m_XScan1 - 2*vec3_t(m_Dx, m_Dy, m_Dz);
      m_X2 = m_XScan2 + 2*vec3_t(m_Dx, m_Dy, m_Dz);
    } else {
      m_X1 *= 0.1;
      m_X2 *= 0.1;
    }
    new_range = (m_X1-m_X2).abs();
  } while (count < 20 && (old_range - new_range)/old_range > 0.2);

  UpdateNodeIndex(m_Grid);
  UpdateCellIndex(m_Grid);
  GuiMainWindow::pointer()->resetXmlDoc();
  GuiMainWindow::pointer()->clearBCs();
  GuiMainWindow::pointer()->addBC(1, BoundaryCondition("imported", "patch"));

  GuiMainWindow::pointer()->resetProgress(" ", 100);
}
