#include "projection_test.h"
#include "surfaceprojection.h"

#include "beziertriangle.h"

#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceAlgorithm()
{
      EG_TYPENAME;
}

void Projection_test::operate()
{
//    project_picked_point();
//   project_all_points();
//   Bezier_test();
//   checkInterpolationGrid();
//    Bezier_circle_test();
  bezierProjectionTest();
}

void Projection_test::project_picked_point()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  vtkIdType id_node = GuiMainWindow::pointer()->getPickedPoint();
  
  int bc_dst = 18;
  
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  
  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());
  vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
  grid->GetPoints()->SetPoint(id_node, x_new.data());
  
  grid->Modified();
}

void Projection_test::project_all_points()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  int bc_src = 18;
  int bc_dst = 18;
  
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals("test");
  
  QVector <vtkIdType> cells;
  QSet <int> bc_src_set;
  bc_src_set.insert(bc_src);
  getSurfaceCells(bc_src_set,cells,grid);
  
  QVector <bool> alreadyprojected(grid->GetNumberOfPoints(),false);
  foreach(vtkIdType id_cell, cells) {
//     qDebug()<<"id_cell="<<id_cell;
    vtkIdType *pts, N_pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for(int i=0;i<N_pts;i++) {
      vtkIdType id_node=pts[i];
      if(!alreadyprojected[id_node]) {
//         qDebug()<<"i="<<i;
        vec3_t x_old;
        grid->GetPoint(id_node, x_old.data());
        vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
        grid->GetPoints()->SetPoint(id_node, x_new.data());
        alreadyprojected[id_node] = true;
      }
    }
  }
  
  grid->Modified();
}

void Projection_test::Bezier_test()
{
//   if (!GuiMainWindow::pointer()->checkSurfProj()) {
//     GuiMainWindow::pointer()->storeSurfaceProjection();
//   }
//   
//   int bc_dst = 18;
//   GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);

  vec3_t X_200(0,0,0);
  vec3_t X_020(1,0,0);
  vec3_t X_002(cos(deg2rad(60)),sin(deg2rad(60)),0);
//   vec3_t X_002(0,1,0);
  
  vec3_t X_011=0.5*(X_020+X_002)+vec3_t(0.5,0.5,0.5);
  vec3_t X_101=0.5*(X_200+X_002)+vec3_t(-0.5,0.5,0.5);
  vec3_t X_110=0.5*(X_200+X_020)+vec3_t(0,-0.5,0.5);
  
/*  vec3_t X_011=0.5*(X_020+X_002);
  vec3_t X_101=0.5*(X_200+X_002);
  vec3_t X_110=0.5*(X_200+X_020);*/
  
  BezierTriangle B(X_200, X_020, X_002, X_011, X_101, X_110);
  B.writeBezierSurface("bezier.vtu");
}

void Projection_test::checkInterpolationGrid()
{
  setBoundaryCodes(GuiMainWindow::pointer()->getAllBoundaryCodes());
  qDebug()<<"getBoundaryCodes()="<<getBoundaryCodes();
//   prepare();
  
  setAllCells();
//   readSettings();
  readVMD();
  
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
  
  updateNodeInfo(true);
  
  
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  int bc_dst;
  bc_dst = 18;
//   updateNodeInfo(true);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals("test");
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeInterpolationGrid("test");
}

void Projection_test::Bezier_circle_test()
{
  int N=10;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, 6*N_cells, 6*N_points);
  
  vtkIdType offset = 0;
  
  for(int i=0;i<6;i++) {
    vec3_t X_200(0,0,0);
    double alpha = i*60;
    vec3_t X_020(cos(deg2rad(alpha)),sin(deg2rad(alpha)),0);
    vec3_t X_002(cos(deg2rad(alpha+60)),sin(deg2rad(alpha+60)),0);
    
    vec3_t h(cos(deg2rad(alpha+30)),sin(deg2rad(alpha+30)),0);
    
//     vec3_t X_011=0.5*(X_020+X_002)+0.25*h;
//     vec3_t X_101=0.5*(X_200+X_002);
//     vec3_t X_110=0.5*(X_200+X_020);
    
    vec3_t X_011=0.5*(X_020+X_002)+0.25*h;
    vec3_t X_101=X_002;
    vec3_t X_110=X_020;
    
    BezierTriangle B(X_200, X_020, X_002, X_011, X_101, X_110);
    offset += addBezierSurface(&B, bezier, offset, N);
  }

  EG_VTKSP(vtkUnstructuredGridWriter,vtu1);
  vtu1->SetFileName("bezier.vtk");
  vtu1->SetInput(bezier);
  vtu1->Write();
  
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu2);
  vtu2->SetFileName("bezier.vtu");
  vtu2->SetDataModeToBinary();
//   vtu2->SetDataModeToAscii();
  vtu2->SetInput(bezier);
  vtu2->Write();
}

int idx_func2(int N, int i, int j)
{
  int offset = -i*(i-2*N-1)/2;
  return offset+j;
}

void Projection_test::bezierProjectionTest()
{
  vec3_t X_200(0,0,0);
  vec3_t X_020(1,0,0);
//   vec3_t X_002(cos(deg2rad(60)),sin(deg2rad(60)),0);
  vec3_t X_002(0,1,0);
  
/*  vec3_t X_011=0.5*(X_020+X_002)+vec3_t( 0.5*cos(deg2rad(30)), 0.5*sin(deg2rad(30)), 0.5);
  vec3_t X_101=0.5*(X_200+X_002)+vec3_t(-0.5*cos(deg2rad(30)), 0.5*sin(deg2rad(30)), 0.5);
  vec3_t X_110=0.5*(X_200+X_020)+vec3_t(0, -0.5, 0.5);*/
  
/*  vec3_t X_011=0.5*(X_020+X_002);
  vec3_t X_101=0.5*(X_200+X_002);
  vec3_t X_110=0.5*(X_200+X_020);*/
  
  vec3_t X_011=0.5*(X_020+X_002)+vec3_t(0.5,0.5,0);
  vec3_t X_101=0.5*(X_200+X_002);
  vec3_t X_110=0.5*(X_200+X_020);
  
  BezierTriangle bezier_triangle(X_200, X_020, X_002, X_011, X_101, X_110);
  bezier_triangle.writeBezierSurface("bezier.vtu");
  
  int N=10;
  int N_cells = (N-1)*(N-1);
  int N_points = (N*N+N)/2;
  qDebug()<<"N_cells="<<N_cells;
  qDebug()<<"N_points="<<N_points;
  
  EG_VTKSP(vtkUnstructuredGrid,bezier);
  allocateGrid(bezier, N_cells, N_points);
  
  vtkIdType offset = 0;
  vtkIdType node_count = 0;
  
  vec3_t origin = bezier_triangle.m_X_200;
  vec3_t ex = bezier_triangle.m_X_020 - bezier_triangle.m_X_200;
  vec3_t ey = bezier_triangle.m_X_002 - bezier_triangle.m_X_200;
  
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_current, bezier, "node_meshdensity_current" );
  
  vtkDoubleArray *vectors = vtkDoubleArray::New();
  vectors->SetName("normals");
  vectors->SetNumberOfComponents(3);
  vectors->SetNumberOfTuples(bezier->GetNumberOfPoints());
  
  for(int i=0;i<N;i++) {
    for(int j=0;j<N-i;j++) {
      double x = i/(double)(N-1);
      double y = j/(double)(N-1);
      vec3_t g_M = origin + x*ex + y*ey;// + vec3_t(0,0,1) + vec3_t(0.5,0,0);
//       vec3_t g_P = bezier_triangle.projectOnQuadraticBezierTriangle(g_M);
//       vec3_t g_P = bezier_triangle.QuadraticBezierTriangle_g(g_M);
      
      vec3_t l_M = bezier_triangle.globalToLocal(g_M);
      vec2_t t_M = vec2_t(l_M[0],l_M[1]);
      vec2_t t_diff = bezier_triangle.fixedPointFunction(t_M,t_M[0],t_M[1]);
      vec3_t l_diff = vec3_t(t_diff[0], t_diff[1], 0);
      vec3_t g_diff = bezier_triangle.localToGlobal(l_diff) - bezier_triangle.m_X_200;
      vtkIdType id_node = offset + node_count;
      node_meshdensity_current->SetValue(id_node, g_diff.abs());
      double n[3];
      n[0]=g_diff[0];
      n[1]=g_diff[1];
      n[2]=g_diff[2];
      vectors->InsertTuple(id_node,n);
      bezier->GetPoints()->SetPoint(id_node, g_M.data());node_count++;
    }
  }
  
  bezier->GetPointData()->SetVectors(vectors);
  vectors->Delete();
  
  int cell_count = 0;
  for(int i=0;i<N-1;i++) {
    for(int j=0;j<N-1-i;j++) {
      vtkIdType pts_triangle1[3];
      pts_triangle1[0]=offset + idx_func2(N, i  ,j  );
      pts_triangle1[1]=offset + idx_func2(N, i+1,j  );
      pts_triangle1[2]=offset + idx_func2(N, i  ,j+1);
      bezier->InsertNextCell(VTK_TRIANGLE,3,pts_triangle1);cell_count++;
      
      if(i+j<N-2) {
        vtkIdType pts_triangle2[3];
        pts_triangle2[0]=offset + idx_func2(N, i+1,j  );
        pts_triangle2[1]=offset + idx_func2(N, i+1,j+1);
        pts_triangle2[2]=offset + idx_func2(N, i  ,j+1);
        bezier->InsertNextCell(VTK_TRIANGLE,3,pts_triangle2);cell_count++;
      }
    }
  }
  
  offset = node_count;
  
  qDebug()<<"node_count="<<node_count;
  qDebug()<<"cell_count="<<cell_count;
  
  EG_VTKSP(vtkXMLUnstructuredGridWriter,vtu2);
  vtu2->SetFileName("bezierProjectionTest.vtu");
  vtu2->SetDataModeToBinary();
//   vtu2->SetDataModeToAscii();
  vtu2->SetInput(bezier);
  vtu2->Write();
}
