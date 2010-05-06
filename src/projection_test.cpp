//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "projection_test.h"
#include "surfaceprojection.h"

#include "beziertriangle.h"
#include "laplacesmoother.h"

#include "vtkUnstructuredGridWriter.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceAlgorithm() {
  EG_TYPENAME;
}

void Projection_test::operate() {
//    project_picked_point();
//   project_all_points();
//   project_all_points2();
//   Bezier_test();
//   checkInterpolationGrid();
//   Bezier_circle_test();
//   bezierFunctionTest();
//   bezierProjectionTest();
//   bezierQuads();
  
  BezierTriangle bezier_triangle;
  for(int i=0; i<11; i++) {
    QString name;
    
    name = "bezier_equi_"+QString::number(i)+"_";
    qWarning()<<name;
    bezier_triangle = specialTriangle(true, i);
    bezier_triangle.setupTriangle();
    bezier_triangle.saveBezierTriangle(name);
    bezierProjectionTest2(bezier_triangle, name);
    qDebug()<<"name="<<name<<" VALID = "<<bezier_triangle.checkControlPoints();
    
    name="bezier_notequi_"+QString::number(i)+"_";
    qWarning()<<name;
    bezier_triangle = specialTriangle(false, i);
    bezier_triangle.setupTriangle();
    bezier_triangle.saveBezierTriangle(name);
    bezierProjectionTest2(bezier_triangle, name);
    qDebug()<<"name="<<name<<" VALID = "<<bezier_triangle.checkControlPoints();
  }

/*  BezierTriangle bezier_triangle;
  bezier_triangle = specialTriangle(false, 0);
  bezier_triangle.setupTriangle();
  bezierProjectionTest2(bezier_triangle, "extrapolation_");*/
}

BezierTriangle Projection_test::specialTriangle(bool equi, int type) {
  vec3_t X_200, X_020, X_002;
  vec3_t X_011, X_101, X_110;

  if (equi) {
    X_200 = vec3_t(0, 0, 0);
    X_020 = vec3_t(1, 0, 0);
    X_002 = vec3_t(cos(deg2rad(60)), sin(deg2rad(60)), 0);
  } else {
    X_200 = vec3_t(0, 0, 0);
    X_020 = vec3_t(1, 0, 0);
    X_002 = vec3_t(0, 1, 0);
  }

  if (type == 0) {
    X_011 = 0.5 * (X_020 + X_002) + vec3_t(0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0.5);
    X_101 = 0.5 * (X_200 + X_002) + vec3_t(-0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0.5);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0.5);
  } else if (type == 1) {
    X_011 = 0.5 * (X_020 + X_002) + vec3_t(0.5, 0.5, 0.5);
    X_101 = 0.5 * (X_200 + X_002) + vec3_t(-0.5, 0, 0.5);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0.5);
  } else if (type == 2) {
    X_011 = 0.5 * (X_020 + X_002);
    X_101 = 0.5 * (X_200 + X_002);
    X_110 = 0.5 * (X_200 + X_020);
  } else if (type == 3) {
    X_011 = 0.5 * (X_020 + X_002) + vec3_t(0.5, 0.5, 0);
    X_101 = 0.5 * (X_200 + X_002);
    X_110 = 0.5 * (X_200 + X_020);
  } else if (type == 4) {
    X_011 = 1./3.*(X_200 + X_020 + X_002);
    X_101 = 1./3.*(X_200 + X_020 + X_002);
    X_110 = 1./3.*(X_200 + X_020 + X_002);
  } else if (type == 5) {
    X_011 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0);
    X_101 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0);
  } else if (type == 6) {
    X_011 = 0.5 * (X_020 + X_002) + vec3_t(0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0);
    X_101 = 0.5 * (X_200 + X_002) + vec3_t(-0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0);
  } else if (type == 7) {
    X_011 = 0.5 * (X_200 + X_020) + vec3_t(0, -1, 1);
    X_101 = 0.5 * (X_200 + X_020) + vec3_t(0, -1, 1);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0);
  } else if (type == 8) {// bad bezier surface!
    X_011 = 0.5 * (X_200 + X_002) + vec3_t(-0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0.5);
    X_101 = 0.5 * (X_020 + X_002) + vec3_t(0.5 * cos(deg2rad(30)), 0.5 * sin(deg2rad(30)), 0.5);
    X_110 = 0.5 * (X_200 + X_020) + vec3_t(0, -0.5, 0.5);
  } else if (type == 9) {// bad bezier surface!
    X_011 = X_200;
    X_101 = X_020;
    X_110 = X_002;
  } else if (type == 10) {// bad bezier surface!
    X_011 = X_200 + vec3_t(0,0,1);
    X_101 = X_020 + vec3_t(0,0,2);
    X_110 = X_002 + vec3_t(0,0,3);
  }
  return BezierTriangle(X_200, X_020, X_002, X_011, X_101, X_110);
}

void Projection_test::Bezier_test() {
  BezierTriangle B0 = specialTriangle(true, 0);
  B0.writeBezierSurface("bezier0", 10);
  BezierTriangle B7 = specialTriangle(true, 7);
  B7.writeBezierSurface("bezier7", 10);
  BezierTriangle B8 = specialTriangle(true, 8);
  B8.writeBezierSurface("bezier8", 10);
}

void Projection_test::project_picked_point() {
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }

  vtkIdType id_node = GuiMainWindow::pointer()->getPickedPoint();

  int bc_dst = 18;

  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(m_Grid);

  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
  m_Grid->GetPoints()->SetPoint(id_node, x_new.data());

  m_Grid->Modified();
}

void Projection_test::project_all_points() {
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }

  int bc_src = 42;
  int bc_dst = 18;

  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(m_Grid);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals("test");

  QVector <vtkIdType> cells;
  QSet <int> bc_src_set;
  bc_src_set.insert(bc_src);
  getSurfaceCells(bc_src_set, cells, m_Grid);

  QVector <bool> alreadyprojected(m_Grid->GetNumberOfPoints(), false);
  foreach(vtkIdType id_cell, cells) {
//     qDebug()<<"id_cell="<<id_cell;
    vtkIdType *pts, N_pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; i++) {
      vtkIdType id_node = pts[i];
      if (!alreadyprojected[id_node]) {
//         qDebug()<<"i="<<i;
        vec3_t x_old;
        m_Grid->GetPoint(id_node, x_old.data());
        vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
        m_Grid->GetPoints()->SetPoint(id_node, x_new.data());
        alreadyprojected[id_node] = true;
      }
    }
  }

  m_Grid->Modified();
}

void Projection_test::checkInterpolationGrid() {
  setBoundaryCodes(GuiMainWindow::pointer()->getAllBoundaryCodes());
  qDebug() << "getBoundaryCodes()=" << getBoundaryCodes();
//   prepare();

  setAllCells();
//   readSettings();
  readVMD();

  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");//node type

  updateNodeInfo(true);


  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }

  int bc_dst;
  bc_dst = 18;
//   updateNodeInfo(true);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(m_Grid);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals("test");
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeInterpolationGrid("test");
}

void Projection_test::Bezier_circle_test() {
  int N = 10;

  MeshPartition new_grid_partition;
  EG_VTKSP(vtkUnstructuredGrid, bezier_first);
  bool first = true;
  
  for (int i = 0; i < 6; i++) {
    vec3_t X_200(0, 0, 0);
    double alpha = i * 60;
    vec3_t X_020(cos(deg2rad(alpha)), sin(deg2rad(alpha)), 0);
    vec3_t X_002(cos(deg2rad(alpha + 60)), sin(deg2rad(alpha + 60)), 0);

    vec3_t h(cos(deg2rad(alpha + 30)), sin(deg2rad(alpha + 30)), 0);

//     vec3_t X_011=0.5*(X_020+X_002)+0.25*h;
//     vec3_t X_101=0.5*(X_200+X_002);
//     vec3_t X_110=0.5*(X_200+X_020);

    vec3_t X_011 = 0.5 * (X_020 + X_002) + 0.25 * h;
    vec3_t X_101 = X_002;
    vec3_t X_110 = X_020;

    // add the local grid
    BezierTriangle B(X_200, X_020, X_002, X_011, X_101, X_110);
    if(first) {
      first = false;
      B.getBezierSurface(bezier_first, N);
      new_grid_partition.setGrid(bezier_first);
      new_grid_partition.setAllCells();
    }
    else {
      EG_VTKSP(vtkUnstructuredGrid, bezier);
      B.getBezierSurface(bezier, N);
      MeshPartition grid_partition(bezier, true);
      new_grid_partition.addPartition(grid_partition);
    }
    
  }

  saveGrid(new_grid_partition.getGrid(), "bezier_all");
}

void Projection_test::bezierFunctionTest() {
  int N = 10;

  BezierTriangle bezier_triangle = specialTriangle(false, 1);
  bezier_triangle.writeBezierSurface("bezier", N);

  int N_cells = (N - 1) * (N - 1);
  int N_points = (N * N + N) / 2;
  qDebug() << "N_cells=" << N_cells;
  qDebug() << "N_points=" << N_points;

  vec2_t toto = vec2_t(0.5, 0.5);
  qDebug() << toto << "->" << bezier_triangle.fixedPointFunction(toto, toto[0], toto[1]);

  EG_VTKSP(vtkUnstructuredGrid, bezier);
  allocateGrid(bezier, N_cells, N_points);

  vtkIdType offset = 0;
  vtkIdType node_count = 0;

  vec3_t origin = bezier_triangle.m_X_200;
  vec3_t ex = bezier_triangle.m_X_020 - bezier_triangle.m_X_200;
  vec3_t ey = bezier_triangle.m_X_002 - bezier_triangle.m_X_200;

  EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, bezier, "node_meshdensity_current");

  vtkDoubleArray *vectors1 = vtkDoubleArray::New();
  vectors1->SetName("diffs");
  vectors1->SetNumberOfComponents(3);
  vectors1->SetNumberOfTuples(bezier->GetNumberOfPoints());

  vtkDoubleArray *vectors2 = vtkDoubleArray::New();
  vectors2->SetName("jacobi");
  vectors2->SetNumberOfComponents(3);
  vectors2->SetNumberOfTuples(bezier->GetNumberOfPoints());

  vtkDoubleArray *vectors3_n = vtkDoubleArray::New();
  vectors3_n->SetName("normals");
  vectors3_n->SetNumberOfComponents(3);
  vectors3_n->SetNumberOfTuples(bezier->GetNumberOfPoints());

  vtkDoubleArray *vectors3_u1 = vtkDoubleArray::New();
  vectors3_u1->SetName("u1");
  vectors3_u1->SetNumberOfComponents(3);
  vectors3_u1->SetNumberOfTuples(bezier->GetNumberOfPoints());

  vtkDoubleArray *vectors3_u2 = vtkDoubleArray::New();
  vectors3_u2->SetName("u2");
  vectors3_u2->SetNumberOfComponents(3);
  vectors3_u2->SetNumberOfTuples(bezier->GetNumberOfPoints());

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N - i; j++) {

      // calculate original mesh point
      double x = i / (double)(N - 1);
      double y = j / (double)(N - 1);
      vec3_t g_M = origin + x * ex + y * ey;// + vec3_t(0,0,1) + vec3_t(0.5,0,0);
      vec3_t g_P = bezier_triangle.projectOnQuadraticBezierTriangle(g_M);
//       vec3_t g_P = bezier_triangle.quadraticBezierTriangle_g(g_M);
//       qDebug()<<"g_M="<<g_M;
      vec2_t t_M = bezier_triangle.global3DToLocal2D(g_M);

      // calculate diff vectors
      vec2_t t_diff = bezier_triangle.fixedPointFunction(t_M, t_M[0], t_M[1]);
      vec3_t g_diff = bezier_triangle.local2DToGlobal3D(t_diff) - bezier_triangle.m_X_200;

//       qDebug()<<"t_diff="<<t_diff;
//       qDebug()<<"g_diff="<<g_diff;

      // calculate tangent vectors
      vec3_t g_center = 1.0 / 3.0 * (bezier_triangle.m_X_200 + bezier_triangle.m_X_020 + bezier_triangle.m_X_002);
      vec2_t t_center = bezier_triangle.global3DToLocal2D(g_center);

      vec2_t displacement = 0.1 * (t_center - t_M);

      vec2_t t_tangent = bezier_triangle.jacobiMatrix(t_M[0], t_M[1]) * displacement;
      vec3_t g_tangent = bezier_triangle.local2DToGlobal3D(t_tangent) - bezier_triangle.m_X_200;

      // calculate normal vectors
      vec3_t g_normal = bezier_triangle.surfaceNormal(t_M, 0);
      vec3_t g_u1 = bezier_triangle.surfaceNormal(t_M, 1);
      vec3_t g_u2 = bezier_triangle.surfaceNormal(t_M, 2);

      // enter the values
      vtkIdType id_node = offset + node_count;
      node_meshdensity_current->SetValue(id_node, g_diff.abs());
      vectors1->InsertTuple(id_node, g_diff.data());
      vectors2->InsertTuple(id_node, g_tangent.data());
      vectors3_n->InsertTuple(id_node, g_normal.data());
      vectors3_u1->InsertTuple(id_node, g_u1.data());
      vectors3_u2->InsertTuple(id_node, g_u2.data());
      bezier->GetPoints()->SetPoint(id_node, g_P.data()); node_count++;
    }
  }

  bezier->GetPointData()->AddArray(vectors1);
  bezier->GetPointData()->AddArray(vectors2);
  bezier->GetPointData()->AddArray(vectors3_n);
  bezier->GetPointData()->AddArray(vectors3_u1);
  bezier->GetPointData()->AddArray(vectors3_u2);
  vectors1->Delete();
  vectors2->Delete();
  vectors3_n->Delete();
  vectors3_u1->Delete();
  vectors3_u2->Delete();

  int cell_count = 0;
  for (int i = 0; i < N - 1; i++) {
    for (int j = 0; j < N - 1 - i; j++) {
      vtkIdType pts_triangle1[3];
      pts_triangle1[0] = offset + trigrid_idx(N, i  , j);
      pts_triangle1[1] = offset + trigrid_idx(N, i + 1, j);
      pts_triangle1[2] = offset + trigrid_idx(N, i  , j + 1);
      bezier->InsertNextCell(VTK_TRIANGLE, 3, pts_triangle1); cell_count++;

      if (i + j < N - 2) {
        vtkIdType pts_triangle2[3];
        pts_triangle2[0] = offset + trigrid_idx(N, i + 1, j);
        pts_triangle2[1] = offset + trigrid_idx(N, i + 1, j + 1);
        pts_triangle2[2] = offset + trigrid_idx(N, i  , j + 1);
        bezier->InsertNextCell(VTK_TRIANGLE, 3, pts_triangle2); cell_count++;
      }
    }
  }

  offset = node_count;

  qDebug() << "node_count=" << node_count;
  qDebug() << "cell_count=" << cell_count;

  saveGrid(bezier, "bezierFunctionTest");
}

void Projection_test::bezierProjectionTest() {
  int N = 10;

  BezierTriangle bezier_triangle = specialTriangle(true, 0);
  bezier_triangle.writeBezierSurface("bezier", N);

  int N_cells = (N - 1) * (N - 1);
  int N_points = (N * N + N) / 2;
  qDebug() << "N_cells=" << N_cells;
  qDebug() << "N_points=" << N_points;

  vec3_t g_center = 1.0 / 3.0 * (bezier_triangle.m_X_200 + bezier_triangle.m_X_020 + bezier_triangle.m_X_002);
  vec2_t t_center = bezier_triangle.global3DToLocal2D(g_center);
  qDebug() << "g_center=" << g_center;
  qDebug() << "t_center=" << t_center;

  vec3_t g_toto = vec3_t(0.5, 0.1, 0);
//   qDebug()<<toto<<"->"<<bezier_triangle.fixedPointFunction(toto,toto[0],toto[1]);
  qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++";
  bezier_triangle.projectOnQuadraticBezierTriangle(g_toto);
  qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++";
//   return;

  EG_VTKSP(vtkUnstructuredGrid, bezier);
  allocateGrid(bezier, N_cells, N_points);

  vtkIdType offset = 0;
  vtkIdType node_count = 0;

  vec3_t origin = bezier_triangle.m_X_200;
  vec3_t ex = bezier_triangle.m_X_020 - bezier_triangle.m_X_200;
  vec3_t ey = bezier_triangle.m_X_002 - bezier_triangle.m_X_200;

  EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, bezier, "node_meshdensity_current");

  vtkDoubleArray *vectors1 = vtkDoubleArray::New();
  vectors1->SetName("normals");
  vectors1->SetNumberOfComponents(3);
  vectors1->SetNumberOfTuples(bezier->GetNumberOfPoints());

  vtkDoubleArray *vectors2 = vtkDoubleArray::New();
  vectors2->SetName("jacobi");
  vectors2->SetNumberOfComponents(3);
  vectors2->SetNumberOfTuples(bezier->GetNumberOfPoints());


  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N - i; j++) {

      // calculate original mesh point
      double x = i / (double)(N - 1);
      double y = j / (double)(N - 1);
      vec3_t g_M = origin + x * ex + y * ey;// + vec3_t(0,0,1) + vec3_t(0.5,0,0);
      vec3_t g_P = bezier_triangle.projectOnQuadraticBezierTriangle(g_M);
//       vec3_t g_P = bezier_triangle.quadraticBezierTriangle_g(g_M);
      qDebug() << "g_M=" << g_M;
      vec2_t t_M = bezier_triangle.global3DToLocal2D(g_M);

      // calculate diff vectors
      vec2_t t_diff = bezier_triangle.fixedPointFunction(t_M, t_M[0], t_M[1]);
      vec3_t g_diff = bezier_triangle.local2DToGlobal3D(t_diff) - bezier_triangle.m_X_200;

      qDebug() << "t_diff=" << t_diff;
      qDebug() << "g_diff=" << g_diff;

      // calculate "jacobi vectors"
      vec2_t displacement = 0.1 * (t_center - t_M);

      vec2_t t_tangent = bezier_triangle.jacobiMatrix(t_M[0], t_M[1]) * displacement;
//       vec2_t t_tangent = bezier_triangle.jacobiMatrix_numeric(t_M,t_M[0],t_M[1],displacement[0],displacement[1]) * displacement;

      vec3_t g_tangent = bezier_triangle.local2DToGlobal3D(t_tangent) - bezier_triangle.m_X_200;

      // enter the values
      vtkIdType id_node = offset + node_count;
      node_meshdensity_current->SetValue(id_node, g_diff.abs());
      vectors1->InsertTuple(id_node, g_diff.data());
      vectors2->InsertTuple(id_node, g_tangent.data());
      bezier->GetPoints()->SetPoint(id_node, g_P.data()); node_count++;
    }
  }

  bezier->GetPointData()->SetVectors(vectors2);
  vectors1->Delete();
  vectors2->Delete();

  int cell_count = 0;
  for (int i = 0; i < N - 1; i++) {
    for (int j = 0; j < N - 1 - i; j++) {
      vtkIdType pts_triangle1[3];
      pts_triangle1[0] = offset + trigrid_idx(N, i  , j);
      pts_triangle1[1] = offset + trigrid_idx(N, i + 1, j);
      pts_triangle1[2] = offset + trigrid_idx(N, i  , j + 1);
      bezier->InsertNextCell(VTK_TRIANGLE, 3, pts_triangle1); cell_count++;

      if (i + j < N - 2) {
        vtkIdType pts_triangle2[3];
        pts_triangle2[0] = offset + trigrid_idx(N, i + 1, j);
        pts_triangle2[1] = offset + trigrid_idx(N, i + 1, j + 1);
        pts_triangle2[2] = offset + trigrid_idx(N, i  , j + 1);
        bezier->InsertNextCell(VTK_TRIANGLE, 3, pts_triangle2); cell_count++;
      }
    }
  }

  offset = node_count;

  qDebug() << "node_count=" << node_count;
  qDebug() << "cell_count=" << cell_count;

  saveGrid(bezier, "bezierProjectionTest");
}

void Projection_test::bezierQuads() {
  int N = 10;

  BezierTriangle bezier_triangle = specialTriangle(true, 0);
  bezier_triangle.writeBezierSurface("bezier", N);

  int N_cells = (N - 1) * (N - 1);
  int N_points = N * N;
  qDebug() << "N_cells=" << N_cells;
  qDebug() << "N_points=" << N_points;

  EG_VTKSP(vtkUnstructuredGrid, bezier);
  allocateGrid(bezier, N_cells, N_points);

  vtkIdType offset = 0;
  vtkIdType node_count = 0;

  vec3_t origin = bezier_triangle.m_X_200;
  vec3_t ex = bezier_triangle.m_X_020 - bezier_triangle.m_X_200;
  vec3_t ey = bezier_triangle.m_X_002 - bezier_triangle.m_X_200;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // calculate original mesh point
      double y = j / (double)(N - 1);
      double x = (i / (double)(N - 1)) * (1.0 - y);
//       double x = (i/(double)(N-1));
      vec3_t g_M = origin + x * ex + y * ey;// + vec3_t(0,0,1) + vec3_t(0.5,0,0);
      vec3_t g_P = bezier_triangle.quadraticBezierTriangle_g(g_M);
      // enter the values
      vtkIdType id_node = offset + node_count;
      bezier->GetPoints()->SetPoint(id_node, g_P.data()); node_count++;
    }
  }

  int cell_count = 0;
  for (int i = 0; i < N - 1; i++) {
    for (int j = 0; j < N - 1; j++) {
      vtkIdType pts_quad[4];
      pts_quad[0] = offset + quadgrid_idx(N, i  , j);
      pts_quad[1] = offset + quadgrid_idx(N, i + 1, j);
      pts_quad[2] = offset + quadgrid_idx(N, i + 1, j + 1);
      pts_quad[3] = offset + quadgrid_idx(N, i  , j + 1);
      bezier->InsertNextCell(VTK_QUAD, 4, pts_quad); cell_count++;
    }
  }

  offset = node_count;

  qDebug() << "node_count=" << node_count;
  qDebug() << "cell_count=" << cell_count;

  saveGrid(bezier, "bezierQuadProjectionTest");
}

// #define EGVTKOBJECT_CREATENODEFIELD(FIELD,TYPE,OW) \
// if (!grid->GetPointData()->GetArray(FIELD)) { \
// EG_VTKSP(TYPE, var); \
// var->SetName(FIELD); \
// var->SetNumberOfValues(Nnodes); \
// grid->GetPointData()->AddArray(var); \
// for (int i = 0; i < grid->GetNumberOfPoints(); ++i) { \
// var->SetValue(i,0); \
// } \
// } else if (OW) { \
// EG_VTKDCN(TYPE, var, grid, FIELD); \
// var->SetNumberOfValues(Nnodes); \
// for (int i = 0; i < grid->GetNumberOfPoints(); ++i) { \
// var->SetValue(i,0); \
// } \
// }

void Projection_test::bezierProjectionTest2(BezierTriangle bezier_triangle, QString prefix) {
  int N = 30;
  bezier_triangle.writeBezierSurface(prefix + "bezier", N);

  bezier_triangle.m_has_neighbour[0] = false;
  bezier_triangle.m_has_neighbour[1] = false;
  bezier_triangle.m_has_neighbour[2] = false;
  bezier_triangle.m_has_neighbour[3] = false;
  bezier_triangle.m_has_neighbour[4] = false;
  bezier_triangle.m_has_neighbour[5] = false;

  int N_cells = (N - 1) * (N - 1);
  int N_points = N * N;
  qDebug() << "N_cells=" << N_cells;
  qDebug() << "N_points=" << N_points;

  EG_VTKSP(vtkUnstructuredGrid, bezier);
  allocateGrid(bezier, N_cells, N_points);

  EG_VTKSP(vtkUnstructuredGrid, bezier_projection);
  allocateGrid(bezier_projection, N_cells, N_points);

  vtkIdType offset = 0;
  vtkIdType node_count = 0;

  vec3_t origin = bezier_triangle.m_X_200;
  vec3_t ex = bezier_triangle.m_X_020 - bezier_triangle.m_X_200;
  vec3_t ey = bezier_triangle.m_X_002 - bezier_triangle.m_X_200;

  int I = 8, J = 6;
  double X = I / (double)(N - 1);
  double Y = J / (double)(N - 1);
  vec3_t g_toto = origin + X * ex + Y * ey;
  qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++";
  qDebug() << "g_toto=" << g_toto << "->" << bezier_triangle.projectOnQuadraticBezierTriangle(g_toto);
  qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++";
//   return;

  vtkDoubleArray *vectors_normals = vtkDoubleArray::New();
  vectors_normals->SetName("normals");
  vectors_normals->SetNumberOfComponents(3);
  vectors_normals->SetNumberOfTuples(bezier_projection->GetNumberOfPoints());

  EG_VTKSP(vtkIntArray, inside0);
  inside0->SetName("inside0");
  inside0->SetNumberOfValues(bezier_projection->GetNumberOfPoints());
  bezier_projection->GetPointData()->AddArray(inside0);
  
  EG_VTKSP(vtkIntArray, inside1);
  inside1->SetName("inside1");
  inside1->SetNumberOfValues(bezier_projection->GetNumberOfPoints());
  bezier_projection->GetPointData()->AddArray(inside1);
  
  EG_VTKSP(vtkIntArray, inside2);
  inside2->SetName("inside2");
  inside2->SetNumberOfValues(bezier_projection->GetNumberOfPoints());
  bezier_projection->GetPointData()->AddArray(inside2);
  
  EG_VTKSP(vtkIntArray, inside_all);
  inside_all->SetName("inside_all");
  inside_all->SetNumberOfValues(bezier_projection->GetNumberOfPoints());
  bezier_projection->GetPointData()->AddArray(inside_all);
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // calculate original mesh point
      double y = -1 + 3 * j / (double)(N - 1);
      double x = -1 + 3 * i / (double)(N - 1);

      /*      double y = j/(double)(N-1);
            double x = i/(double)(N-1);*/

//       double x = (i/(double)(N-1));
      vec3_t g_M = origin + x * ex + y * ey;// + vec3_t(0,0,1) + vec3_t(0.5,0,0);
      vec2_t t_M = bezier_triangle.global3DToLocal2D(g_M);

      vec3_t g_P = bezier_triangle.quadraticBezierTriangle_g(g_M);
      vec3_t g_P_projection = bezier_triangle.projectOnQuadraticBezierTriangle(g_M, 0);
      double L,u;
//       vec3_t g_P_projection = bezier_triangle.projectOnBezierSide(g_M,2,L,u);
      
//       vec3_t g_normal = bezier_triangle.surfaceNormal(t_M, 0);
      vec3_t g_normal = bezier_triangle.projectOnQuadraticBezierTriangle(g_M, 1);
//       vec3_t g_normal(0,0,0);
      vec2_t t_tangent;
      bool I0 = bezier_triangle.insideBezierCurve(t_M,0,t_tangent);
      bool I1 = bezier_triangle.insideBezierCurve(t_M,1,t_tangent);
      bool I2 = bezier_triangle.insideBezierCurve(t_M,2,t_tangent);
      bool I_all = bezier_triangle.insideBezierSurface(g_M);
      
      // enter the values
      vtkIdType id_node = offset + node_count;

      bezier->GetPoints()->SetPoint(id_node, g_P.data());
      bezier_projection->GetPoints()->SetPoint(id_node, g_P_projection.data());
      vectors_normals->InsertTuple(id_node, g_normal.data());
      inside0->SetValue(id_node,I0);
      inside1->SetValue(id_node,I1);
      inside2->SetValue(id_node,I2);
      inside_all->SetValue(id_node,I_all);
      node_count++;
    }
  }

  bezier_projection->GetPointData()->AddArray(vectors_normals);
  vectors_normals->Delete();

  int cell_count = 0;
  for (int i = 0; i < N - 1; i++) {
    for (int j = 0; j < N - 1; j++) {
      vtkIdType pts_quad[4];
      pts_quad[0] = offset + quadgrid_idx(N, i  , j);
      pts_quad[1] = offset + quadgrid_idx(N, i + 1, j);
      pts_quad[2] = offset + quadgrid_idx(N, i + 1, j + 1);
      pts_quad[3] = offset + quadgrid_idx(N, i  , j + 1);

      bezier->InsertNextCell(VTK_QUAD, 4, pts_quad);
      bezier_projection->InsertNextCell(VTK_QUAD, 4, pts_quad);
      cell_count++;
    }
  }

  offset = node_count;

  qDebug() << "node_count=" << node_count;
  qDebug() << "cell_count=" << cell_count;

  saveGrid(bezier, prefix + "bezierQuadTest");
  saveGrid(bezier_projection, prefix + "QuadProjectionTest");
}

void Projection_test::project_all_points2() {
  
  // store projection surface
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  writeGrid(m_Grid, "allinone_before_meshing");
  
  // mesh surface
  SurfaceMesher mesher;
  mesher();
  
  writeGrid(m_Grid, "allinone_after_meshing");
  
  // correct curvature
  LaplaceSmoother lap;
  lap.setCorrectCurvature(true);
  lap.setNoCheck(true);
  lap.setGrid(m_Grid);
  QVector<vtkIdType> cls;
  getSurfaceCells(m_BoundaryCodes, cls, m_Grid);
  lap.setCells(cls);
  lap.setNumberOfIterations(1);
  m_BoundaryCodes = GuiMainWindow::pointer()->getAllBoundaryCodes();
  lap.setBoundaryCodes(m_BoundaryCodes);//IMPORTANT: so that unselected nodes become fixed when node types are updated!
  if (m_UseProjectionForSmoothing) {
    lap.setProjectionOn();
  } else {
    lap.setProjectionOff();
  }
  if (m_UseNormalCorrectionForSmoothing) {
    lap.setNormalCorrectionOn();
  } else {
    lap.setNormalCorrectionOff();
  }
  lap();
  m_SmoothSuccess = lap.succeeded();
  
  writeGrid(m_Grid, "allinone_after_correctCurvature");
  
  return;
  EG_BUG;
  
  QSet<int> bcs = GuiMainWindow::pointer()->getAllBoundaryCodes();
  foreach (int bc, bcs) {
    GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(m_Grid);
    GuiMainWindow::pointer()->getSurfProj(bc)->setCorrectCurvature(true);
  }
  
  QVector <bool> alreadyprojected(m_Grid->GetNumberOfPoints(), false);
  QVector <vtkIdType> cells;
  getAllSurfaceCells(cells, m_Grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  
  foreach(vtkIdType id_cell, cells) {
    int bc = cell_code->GetValue(id_cell);
    
    vtkIdType *pts, N_pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    for (int i = 0; i < N_pts; i++) {
      vtkIdType id_node = pts[i];
      if (!alreadyprojected[id_node]) {
        vec3_t x_old;
        m_Grid->GetPoint(id_node, x_old.data());
        vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_old, id_node);
        m_Grid->GetPoints()->SetPoint(id_node, x_new.data());
        alreadyprojected[id_node] = true;
      }
    }
  }
  
  m_Grid->Modified();
}

bool Projection_test::setNewPosition(vtkIdType id_node, vec3_t x_new)
{
  using namespace GeometryTools;
  
  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  m_Grid->GetPoints()->SetPoint(id_node, x_new.data());
  bool move = true;
  
  vec3_t n(0,0,0);
  QVector<vec3_t> cell_normals(m_Part.n2cGSize(id_node));
  double A_max = 0;//area of the biggest neighbour cell of id_node
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    double A = fabs(GeometryTools::cellVA(m_Grid, m_Part.n2cGG(id_node, i)));
    A_max = max(A, A_max);
    cell_normals[i] = GeometryTools::cellNormal(m_Grid, m_Part.n2cGG(id_node, i));
    cell_normals[i].normalise();
  }
  int N = 0;
  for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
    double A = fabs(GeometryTools::cellVA(m_Grid, m_Part.n2cGG(id_node, i)));
    if (A > 0.01*A_max) {
      n += cell_normals[i];
      ++N;
    }
  }
  if (N == 0) {
    EG_BUG;
    move = false;
  } else {
    n.normalise();
    double L_max = 0;// distance to the furthest neighbour node of id_node
    for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
      vec3_t xn;
      m_Grid->GetPoint(m_Part.n2nGG(id_node, i), xn.data());
      double L = (xn - x_old).abs();
      L_max = max(L, L_max);
    }
    vec3_t x_summit = x_new + L_max*n;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vec3_t x[3];
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(m_Part.n2cGG(id_node, i), N_pts, pts);
      if (N_pts != 3) {
        EG_BUG;
      }
      for (int j = 0; j < N_pts; ++j) {
        m_Grid->GetPoint(pts[j], x[j].data());
      }
      if (GeometryTools::tetraVol(x[0], x[1], x[2], x_summit, false) <= 0) {
//         saveGrid(m_Grid, "after_move1");
//         qWarning()<<"Cannot move point "<<id_node<<" because negative tetraVol.";
        move = false;
//         m_Grid->GetPoints()->SetPoint(id_node, x_old.data());
//         saveGrid(m_Grid, "before_move1");
//         EG_BUG;
        break;
      }
    }
  }
  if (move) {
    for (int i = 0; i < cell_normals.size(); ++i) {
      for (int j = 0; j < cell_normals.size(); ++j) {
        if (cell_normals[i]*cell_normals[j] < -1000*0.1) {//Why -1000*0.1?
          saveGrid(m_Grid, "after_move2");
          qWarning()<<"Cannot move point "<<id_node<<" because opposite cell_normals.";
          move = false;
          m_Grid->GetPoints()->SetPoint(id_node, x_old.data());
          saveGrid(m_Grid, "before_move2");
          EG_BUG;
          break;
        }
      }
    }
  }
  
  if (!move) {
    // comment this out if you want points to always move
    m_Grid->GetPoints()->SetPoint(id_node, x_old.data());
//     saveGrid(m_Grid, "before_move");
  }
  return move;
}

bool Projection_test::moveNode(vtkIdType id_node, vec3_t &Dx)
{
  return false;
/*  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  bool moved = false;
  for (int i_relaxation = 0; i_relaxation < 1; ++i_relaxation) {
    vec3_t x_new = x_old + Dx;
    int i_nodes = m_Part.localNode(id_node);
    if (m_NodeToBc[i_nodes].size() == 1) {
      int bc = m_NodeToBc[i_nodes][0];
      x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new, id_node);
    } else {
      for (int i_proj_iter = 0; i_proj_iter < 20; ++i_proj_iter) {
        foreach (int bc, m_NodeToBc[i_nodes]) {
          x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_new, id_node);
        }
      }
    }
    if (setNewPosition(id_node, x_new)) {
      moved = true;
      Dx = x_new - x_old;
      break;
    }
    Dx *= 0.5;
  }
  return moved;*/
}
