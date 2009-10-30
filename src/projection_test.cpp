#include "projection_test.h"
#include "surfaceprojection.h"

#include "beziertriangle.h"

#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceOperation()
{
      EG_TYPENAME;
}

void Projection_test::operate()
{
//    project_picked_point();
//   project_all_points();
//   Bezier_test();
  checkInterpolationGrid();
//    Bezier_circle_test();
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
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals();
  
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
//   GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals();

  vec3_t X_200(0,0,0);
  vec3_t X_020(1,0,0);
  vec3_t X_002(cos(deg2rad(60)),sin(deg2rad(60)),0);
  vec3_t X_011=0.5*(X_020+X_002)+vec3_t(0.5,0.5,0.5);
  vec3_t X_101=0.5*(X_200+X_002)+vec3_t(-0.5,0.5,0.5);
  vec3_t X_110=0.5*(X_200+X_020)+vec3_t(0,-0.5,0.5);
  
  BezierTriangle B(X_200, X_020, X_002, X_011, X_101, X_110);
  B.writeBezierSurface();
}

void Projection_test::checkInterpolationGrid()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  int bc_dst;
  bc_dst = 18;
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->writeGridWithNormals();
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setupInterpolationGrid();
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
    
    vec3_t X_011=0.5*(X_020+X_002)+0.25*h;
    vec3_t X_101=0.5*(X_200+X_002);
    vec3_t X_110=0.5*(X_200+X_020);
    
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
