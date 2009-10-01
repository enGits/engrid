#include "projection_test.h"
#include "surfaceprojection.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceOperation()
{
      EG_TYPENAME;
}

void Projection_test::operate()
{
  project_picked_point();
//   project_all_points();
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
  
  int bc_src = 42;
  int bc_dst = 18;
  
  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  
  QVector <vtkIdType> cells;
  QSet <int> bc_src_set;
  bc_src_set.insert(bc_src);
  getSurfaceCells(bc_src_set,cells,grid);
  
  foreach(vtkIdType id_cell, cells) {
    vtkIdType *pts, N_pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for(int i=0;i<N_pts;i++) {
      vtkIdType id_node=pts[i];
      vec3_t x_old;
      grid->GetPoint(id_node, x_old.data());
      vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
      grid->GetPoints()->SetPoint(id_node, x_new.data());
    }
  }
  
  grid->Modified();
}
