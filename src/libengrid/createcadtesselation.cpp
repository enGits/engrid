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
  m_NumIterations = 0;
  m_PreservationType = 0;
  m_SmallestFeatureSize = 1e99;
}

void CreateCadTesselation::scan(bool create_grid, int interlaces)
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
        if (preserveFluid()) {
          g->SetValue(getIdx(i,j,k), 1);
        } else {
          g->SetValue(getIdx(i,j,k), 0);
        }
      }
    }
  }
  m_XScan1 = m_X2;
  m_XScan2 = m_X1;

  vec3_t x_in, x_out, n_in, n_out;

  int progress = 0;
  m_GeometryFound = false;

  double dxi = m_Dx/(interlaces + 1);
  double dyi = m_Dy/(interlaces + 1);
  double dzi = m_Dz/(interlaces + 1);

  double shift = 1e-6;

  // scan ij plane -> k direction
  //      xy plane -> z direction
  for (int i = 0; i < m_Ni; ++i) {
    for (int j = 0; j < m_Nj; ++j) {
      vec3_t x0 = getX(i,j,0);
      for (int i_il = 0; i_il <= interlaces; ++i_il) {
        for (int j_il = 0; j_il <= interlaces; ++j_il) {
          vec3_t x = x0;
          x[0] += i_il*dxi;
          x[1] += j_il*dyi;
          int k_last = 0;
          while (shootRay(x, vec3_t(0,0,1), x_in, x_out, n_in, n_out)) {
            m_GeometryFound = true;
            m_XScan1[2] = min(m_XScan1[2], x_in[2]);
            m_XScan2[2] = max(m_XScan2[2], x_out[2]);
            int k1 = int((x_in[2] - m_X1[2])/m_Dz) + 1;
            int k2 = int((x_out[2] - m_X1[2])/m_Dz);
            if (preserveFluid()) {
              for (int k = k_last; k < min(m_Nk-1,k1); ++k) {
                g->SetValue(getIdx(i,j,k), 0);
              }
            } else {
              for (int k = max(0,k1); k <= min(m_Nk-1,k2); ++k) {
                g->SetValue(getIdx(i,j,k), 1);
              }
            }
            x = x_out;
            x[2] += shift*m_Dz;
            k_last = min(m_Nk-1,k2+1);
          }
          if (preserveFluid()) {
            for (int k = k_last; k <= m_Nk-1; ++k) {
              g->SetValue(getIdx(i,j,k), 0);
            }
          }
        }
      }
    }
    progress += m_Nj;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  // scan ik plane -> j direction
  //      xz plane -> y direction
  for (int i = 0; i < m_Ni; ++i) {
    for (int k = 0; k < m_Nk; ++k) {
      vec3_t x0 = getX(i,0,k);
      for (int i_il = 0; i_il <= interlaces; ++i_il) {
        for (int k_il = 0; k_il <= interlaces; ++k_il) {
          vec3_t x = x0;
          x[0] += i_il*dxi;
          x[2] += k_il*dzi;
          int j_last = 0;
          /*
          if (fabs(x0[0]) < 1.5 && fabs(x0[2]) < 1.5 && create_grid) {
            cout << x0 << endl;
          }
          */
          while (shootRay(x, vec3_t(0,1,0), x_in, x_out, n_in, n_out)) {
            m_GeometryFound = true;
            m_XScan1[1] = min(m_XScan1[1], x_in[1]);
            m_XScan2[1] = max(m_XScan2[1], x_out[1]);
            int j1 = int((x_in[1]-m_X1[1])/m_Dy) + 1;
            int j2 = int((x_out[1]-m_X1[1])/m_Dy);
            if (preserveFluid()) {
              for (int j = j_last; j < min(m_Nj-1,j1); ++j) {
                g->SetValue(getIdx(i,j,k), 0);
              }
            } else {
              for (int j = max(0,j1); j <= min(m_Nj-1,j2); ++j) {
                g->SetValue(getIdx(i,j,k), 1);
              }
            }
            x = x_out;
            x[1] += shift*m_Dy;
            j_last = min(m_Nj-1,j2+1);
          }
          if (preserveFluid()) {
            for (int j = j_last; j <= m_Nj-1; ++j) {
              g->SetValue(getIdx(i,j,k), 0);
            }
          }
        }
      }
    }
    progress += m_Nk;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  // scan jk plane -> i direction
  //      yz plane -> x direction
  for (int j = 0; j < m_Nj; ++j) {
    for (int k = 0; k < m_Nk; ++k) {
      vec3_t x0 = getX(0,j,k);
      for (int j_il = 0; j_il <= interlaces; ++j_il) {
        for (int k_il = 0; k_il <= interlaces; ++k_il) {
          vec3_t x = x0;
          x[1] += j_il*dyi;
          x[2] += k_il*dzi;
          int i_last = 0;
          while (shootRay(x, vec3_t(1,0,0), x_in, x_out, n_in, n_out)) {
            m_GeometryFound = true;
            m_XScan1[0] = min(m_XScan1[0], x_in[0]);
            m_XScan2[0] = max(m_XScan2[0], x_out[0]);
            int i1 = int((x_in[0]-m_X1[0])/m_Dx) + 1;
            int i2 = int((x_out[0]-m_X1[0])/m_Dx);
            if (preserveFluid()) {
              for (int i = i_last; i < min(m_Ni-1,i1); ++i) {
                g->SetValue(getIdx(i,j,k), 0);
              }
            } else {
              for (int i = max(0,i1); i <= min(m_Ni-1,i2); ++i) {
                g->SetValue(getIdx(i,j,k), 1);
              }
            }
            x = x_out;
            x[0] += shift*m_Dx;
            i_last = min(m_Ni-1,i2+1);
          }
          if (preserveFluid()) {
            for (int i = i_last; i <= m_Ni-1; ++i) {
              g->SetValue(getIdx(i,j,k), 0);
            }
          }
        }
      }
    }
    progress += m_Nk;
    GuiMainWindow::pointer()->setProgress(progress);
  }

  if (create_grid) {

    GuiMainWindow::pointer()->resetProgress("smoothing", m_Ni*m_Nj*m_Nk*m_NumIterations);

    // smooth image data
    int progress = 0;
    for (int iter = 1; iter <= m_NumIterations; ++iter) {
      for (int i = 1; i < m_Ni-1; ++i) {
        for (int j = 1; j < m_Nj-1; ++j) {
          for (int k = 1; k < m_Nk-1; ++k) {
            int idx = getIdx(i,j,k);
            bool preserve = false;
            if (m_PreservationType == 1 && g->GetValue(idx) > 0.999) {
              preserve = true;
            }
            if (m_PreservationType == 2 && g->GetValue(idx) < 0.001) {
              preserve = true;
            }
            if (!preserve) {
              double g_new = 0;
              g_new += g->GetValue(getIdx(i+1,j,k));
              g_new += g->GetValue(getIdx(i-1,j,k));
              g_new += g->GetValue(getIdx(i,j+1,k));
              g_new += g->GetValue(getIdx(i,j-1,k));
              g_new += g->GetValue(getIdx(i,j,k+1));
              g_new += g->GetValue(getIdx(i,j,k-1));
              g_new *= 1.0/6.0;
              g->SetValue(idx, g_new);
            }
          }
        }
        progress += m_Nj*m_Nk;
        GuiMainWindow::pointer()->setProgress(progress);
      }
    }    

    // write image data for testing and reporting purposes
    EG_VTKSP(vtkXMLImageDataWriter, vti);
    QString file_name = GuiMainWindow::pointer()->getCwd() + "/g.vti";
    vti->SetFileName(qPrintable(file_name));
    vti->SetDataModeToBinary();
    vti->SetInput(gdata);
    vti->Write();

    gdata->GetPointData()->SetActiveScalars("g");
    EG_VTKSP(vtkContourFilter, contour);
    contour->SetInput(gdata);
    contour->SetNumberOfContours(1);
    double g_level = 0.5;
    if (m_PreservationType == 1) {
      g_level = 0.5;
    }
    if (m_PreservationType == 2) {
      g_level = 0.5;
    }
    contour->SetValue(0, g_level);
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
      x[0] = m_X1[0] + x[0]*m_Dx;
      x[1] = m_X1[1] + x[1]*m_Dx;
      x[2] = m_X1[2] + x[2]*m_Dx;
      m_Grid->GetPoints()->SetPoint(id_node, x.data());
    }
    EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
    for (vtkIdType id_cell = 0; id_cell < decimate->GetOutput()->GetNumberOfPolys(); ++id_cell) {
      EG_GET_CELL(id_cell, decimate->GetOutput());
      vtkIdType id_new_cell = m_Grid->InsertNextCell(type_cell, num_pts, pts);
      cell_code->SetValue(id_new_cell, 1);
    }
  }
}

void CreateCadTesselation::operate()
{
  double preview_memory = min(1024*1024*100.0, m_ScanMemory);
  int N = int(pow(preview_memory/sizeof(float), 1.0/3.0));
  m_Ni = N;
  m_Nj = N;
  m_Nk = N;
  m_X1 = vec3_t(-1e6,-1e6,-1e6);
  m_X2 = vec3_t( 1e6, 1e6, 1e6);
  m_Dx = (m_X2[0] - m_X1[0])/(m_Ni-1);
  m_Dy = (m_X2[1] - m_X1[1])/(m_Nj-1);
  m_Dz = (m_X2[2] - m_X1[2])/(m_Nk-1);
  double new_vol, old_vol;
  int count = 0;
  do {
    old_vol = (m_X2[0]-m_X1[0])*(m_X2[1]-m_X1[1])*(m_X2[2]-m_X1[2]);
    QString num;
    QString text = "scanning (V=";
    num.setNum(old_vol);
    text += num + ")";

    GuiMainWindow::pointer()->resetProgress(text, 3*N*N);
    scan(false);
    if (m_GeometryFound) {
      m_X1 = m_XScan1 - 2*vec3_t(m_Dx, m_Dy, m_Dz);
      m_X2 = m_XScan2 + 2*vec3_t(m_Dx, m_Dy, m_Dz);
    } else {
      m_X1 *= 0.1;
      m_X2 *= 0.1;
    }
    new_vol = (m_X2[0]-m_X1[0])*(m_X2[1]-m_X1[1])*(m_X2[2]-m_X1[2]);
  } while (count < 20 && (old_vol - new_vol)/old_vol > 0.05);

  // bounding box should now be established
  // last scan run with the full resoluion (if required)
  double Lx = m_X2[0] - m_X1[0];
  double Ly = m_X2[1] - m_X1[1];
  double Lz = m_X2[2] - m_X1[2];
  double max_size = m_ScanMemory/sizeof(float);
  double delta = pow(Lx*Ly*Lz/max_size, 1.0/3.0);
  int interlaces = 0;
  if (preserveFluid() || preserveSolid()) {
    interlaces = int(2*delta/m_SmallestFeatureSize);
  }
  m_Ni = int(Lx/delta) + 1;
  m_Nj = int(Ly/delta) + 1;
  m_Nk = int(Lz/delta) + 1;
  QString num;
  QString text = "scanning (h=";
  num.setNum(delta/(interlaces+1));
  text += num + ")";
  GuiMainWindow::pointer()->resetProgress(text, m_Ni*m_Nj + m_Ni*m_Nk + m_Nj*m_Nk);
  scan(true, interlaces);

  UpdateNodeIndex(m_Grid);
  UpdateCellIndex(m_Grid);
  GuiMainWindow::pointer()->resetXmlDoc();
  GuiMainWindow::pointer()->clearBCs();
  GuiMainWindow::pointer()->addBC(1, BoundaryCondition("imported", "patch"));

  GuiMainWindow::pointer()->resetProgress(" ", 100);
}
