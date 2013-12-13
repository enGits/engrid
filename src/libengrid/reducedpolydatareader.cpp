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

#include "reducedpolydatareader.h"
#include "guimainwindow.h"
#include "octree.h"

#include <QFileInfo>

#include <vtkPolyDataReader.h>
#include <vtkEgPolyDataToUnstructuredGridFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkClipVolume.h>
#include <vtkContourFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>

ReducedPolyDataReader::ReducedPolyDataReader()
{
  setFormat("legacy VTK files(*.vtk)");
  setExtension(".vtk");
  m_MaxEdgeLength = 10;
}

void ReducedPolyDataReader::computeLevelSet(vtkUnstructuredGrid* m_Grid, vtkPolyData* poly)
{
  // create triangles
  poly->BuildCells();
  QVector<Triangle> triangles(poly->GetNumberOfPolys());
  for (vtkIdType id_poly = 0; id_poly < poly->GetNumberOfPolys(); ++id_poly) {
    vtkIdType Npts, *pts;
    poly->GetCellPoints(id_poly, Npts, pts);
    if (Npts == 3) {
      poly->GetPoint(pts[0], triangles[id_poly].a.data());
      poly->GetPoint(pts[1], triangles[id_poly].b.data());
      poly->GetPoint(pts[2], triangles[id_poly].c.data());
      triangles[id_poly].id_a = pts[0];
      triangles[id_poly].id_b = pts[1];
      triangles[id_poly].id_c = pts[2];
      triangles[id_poly].g1 = triangles[id_poly].b - triangles[id_poly].a;
      triangles[id_poly].g2 = triangles[id_poly].c - triangles[id_poly].a;
      triangles[id_poly].g3 = triangles[id_poly].g1.cross(triangles[id_poly].g2);
      triangles[id_poly].A  = 0.5*triangles[id_poly].g3.abs();
      triangles[id_poly].g3.normalise();
      triangles[id_poly].G.column(0, triangles[id_poly].g1);
      triangles[id_poly].G.column(1, triangles[id_poly].g2);
      triangles[id_poly].G.column(2, triangles[id_poly].g3);
      triangles[id_poly].GI = triangles[id_poly].G.inverse();
      triangles[id_poly].smallest_length = (triangles[id_poly].b - triangles[id_poly].a).abs();
      triangles[id_poly].smallest_length = min(triangles[id_poly].smallest_length, (triangles[id_poly].c - triangles[id_poly].b).abs());
      triangles[id_poly].smallest_length = min(triangles[id_poly].smallest_length, (triangles[id_poly].a - triangles[id_poly].c).abs());
    } else {
      EG_BUG;
    }
  }

  //vtkDoubleArray *scalars = vtkDoubleArray::New();
  EG_VTKSP(vtkDoubleArray, scalars);
  scalars->SetNumberOfComponents(1);
  scalars->SetName("g");
  scalars->Allocate(m_Grid->GetNumberOfPoints());

  QProgressDialog progress("Reducing triangulation", "Abort", 0, m_Grid->GetNumberOfPoints());
  progress.setWindowModality(Qt::ApplicationModal);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    progress.setValue(id_node);
    QApplication::processEvents();
    if (progress.wasCanceled()) {
      EG_ERR_RETURN("interrupted by user");
    }
    double g_levelset = 1e99;
    vec3_t xp;
    m_Grid->GetPoint(id_node, xp.data());
    foreach (Triangle T, triangles) {
      vec3_t xi(1e99,1e99,1e99);
      vec3_t ri;
      double scal = (xp - T.a)*T.g3;
      vec3_t x1, x2;
      if (scal > 0) {
        x1 = xp + T.g3;
        x2 = xp - scal*T.g3 - T.g3;
      } else {
        x1 = xp - T.g3;
        x2 = xp - scal*T.g3 + T.g3;
      }
      double d = 1e99;

      bool intersects_face = GeometryTools::intersectEdgeAndTriangle(T.a, T.b, T.c, x1, x2, xi, ri);
      if (intersects_face) {
        vec3_t dx = xp - xi;
        d = fabs(dx*T.g3);
      } else {
        double kab = GeometryTools::intersection(T.a, T.b - T.a, xp, T.b - T.a);
        double kac = GeometryTools::intersection(T.a, T.c - T.a, xp, T.c - T.a);
        double kbc = GeometryTools::intersection(T.b, T.c - T.b, xp, T.c - T.b);
        double dab = (T.a + kab*(T.b-T.a) - xp).abs();
        double dac = (T.a + kac*(T.c-T.a) - xp).abs();
        double dbc = (T.b + kbc*(T.c-T.b) - xp).abs();
        bool set = false;
        if ((kab >= 0) && (kab <= 1)) {
          if (dab < d) {
            xi = T.a + kab*(T.b-T.a);
            d = dab;
            set = true;
          }
        }
        if ((kac >= 0) && (kac <= 1)) {
          if (dac < d) {
            xi = T.a + kac*(T.c-T.a);
            d = dac;
            set = true;
          }
        }
        if ((kbc >= 0) && (kbc <= 1)) {
          if (dbc < d) {
            xi = T.b + kbc*(T.c-T.b);
            d = dbc;
            set = true;
          }
        }
        double da = (T.a - xp).abs();
        double db = (T.b - xp).abs();
        double dc = (T.c - xp).abs();
        if (da < d) {
          xi = T.a;
          d = da;
          set = true;
        }
        if (db < d) {
          xi = T.b;
          d = db;
        }
        if (dc < d) {
          xi = T.c;
          d = dc;
          set = true;
        }
        if (!set) {
          EG_BUG;
        }
      }
      if (xi[0] > 1e98) {
        EG_BUG;
      }
      vec3_t dx = xp - xi;
      if (d < fabs(g_levelset)) {
        if (dx*T.g3 > 0) {
          g_levelset = d;
        } else {
          g_levelset = -d;
        }
      }
    }
    scalars->InsertTuple1(id_node, g_levelset);
  }
  progress.setValue(m_Grid->GetNumberOfPoints());
  m_Grid->GetPointData()->SetScalars(scalars);
}


void ReducedPolyDataReader::operate()
{
  try {
    QFileInfo file_info(GuiMainWindow::pointer()->getFilename());
    readInputFileName(file_info.completeBaseName() + ".vtk");
    if (isValid()) {
      EG_VTKSP(vtkPolyDataReader, vtk);
      double bounds[6];
      vtk->SetFileName(qPrintable(getFileName()));
      vtk->Update();
      vtk->GetOutput()->GetBounds(bounds);
      vec3_t x1(bounds[0], bounds[2], bounds[4]);
      vec3_t x2(bounds[1], bounds[3], bounds[5]);
      double m_MaxEdgeLength = (x1-x2).abs()/20;
//      double m_MinEdgeLength = (x1-x2).abs()/100;
      vec3_t x12 = 0.5*(x1 + x2);
      vec3_t dx = x2 - x1;
      double d = max(dx[0], max(dx[1], dx[2]));
      dx = vec3_t(d,d,d);
      x1 = x12 - 0.75*dx;
      x2 = x12 + 0.75*dx;
      Octree octree;
      octree.setBounds(x1, x2);
      octree.setSmoothTransitionOn();

      EG_VTKSP(vtkSmoothPolyDataFilter, smooth1);
      smooth1->SetInputConnection(vtk->GetOutputPort());
      smooth1->FeatureEdgeSmoothingOn();
      smooth1->BoundarySmoothingOn();
      smooth1->SetNumberOfIterations(200);
      smooth1->Update();
      smooth1->GetOutput()->BuildCells();

      /*
      QVector<double> max_length(smooth->GetOutput()->GetNumberOfPoints(), m_MaxEdgeLength);
      for (vtkIdType id_node = 0; id_node < smooth->GetOutput()->GetNumberOfPolys(); ++id_node) {
        vtkIdType N_pts, *pts;
        smooth1->GetOutput()->GetCellPoints(id_node, N_pts, pts);
        vec3_t x[N_pts];
        for (int i = 0; i < N_pts; ++i) {
          smooth1->GetOutput()->GetPoint(pts[i], x[i].data());
        }
        for (int i = 0; i < N_pts-1; ++i) {
          int j = i+1;
          if (j >= N_pts) {
            j = 0;
          }
          double L = (x[i]-x[j]).abs();
          max_length[pts[i]] = max(max_length[pts[i]], L);
          max_length[pts[j]] = max(max_length[pts[j]], L);
        }
      }
      */

      EG_VTKSP(vtkDecimatePro, decimate);
      decimate->SetFeatureAngle(20.0);
      decimate->PreserveTopologyOn();
      decimate->SetInputConnection(smooth1->GetOutputPort());
      decimate->SetTargetReduction(0.9);
      EG_VTKSP(vtkSmoothPolyDataFilter, smooth2);
      smooth2->SetInputConnection(decimate->GetOutputPort());
      smooth2->FeatureEdgeSmoothingOn();
      smooth2->BoundarySmoothingOn();
      smooth2->SetNumberOfIterations(200);
      smooth2->Update();

      int N;
      do {
        for (vtkIdType id_node = 0; id_node < smooth2->GetOutput()->GetNumberOfPoints(); ++id_node) {
          vec3_t x;
          smooth2->GetOutput()->GetPoint(id_node, x.data());
          int cell = octree.findCell(x);
          if (octree.getDx(cell) >m_MaxEdgeLength) {
            octree.markToRefine(cell);
          }
        }
        N = octree.refineAll();
      } while (N > 0);

      EG_VTKSP(vtkUnstructuredGrid, oct_grid);
      octree.toVtkGridHangingNodes(oct_grid);

      computeLevelSet(oct_grid, smooth2->GetOutput());

      EG_VTKSP(vtkDelaunay3D, delaunay);
      delaunay->SetInputData(oct_grid);

      EG_VTKSP(vtkContourFilter, contour);
      contour->SetInputConnection(delaunay->GetOutputPort());
      contour->GenerateValues(1, 0, 0);
      contour->SetInputData(oct_grid);
      contour->Update();

      {
        EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu);
        QString file_name = GuiMainWindow::pointer()->getCwd() + "/oct_grid.vtu";
        vtu->SetFileName(qPrintable(file_name));
        vtu->SetDataModeToBinary();
        vtu->SetInputData(oct_grid);
        vtu->Write();
      }
      if (contour->GetOutput()->GetNumberOfPolys() == 0) {
        EG_ERR_RETURN("unable to extract reduced surface");
      }

      EG_VTKSP(vtkEgPolyDataToUnstructuredGridFilter, poly2ug);
      poly2ug->SetInputConnection(contour->GetOutputPort());
      poly2ug->Update();

      makeCopy(poly2ug->GetOutput(), m_Grid);
      createBasicFields(m_Grid, m_Grid->GetNumberOfCells(), m_Grid->GetNumberOfPoints());
      UpdateNodeIndex(m_Grid);
      UpdateCellIndex(m_Grid);
      EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
      for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
        bc->SetValue(id_cell, 0);
      }
    }
  } catch (Error err) {
    err.display();
  }
}

