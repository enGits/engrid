#include "projection_test.h"
#include "surfaceprojection.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceOperation()
{
      EG_TYPENAME;
}

void Projection_test::operate()
{
  if (!GuiMainWindow::pointer()->checkSurfProj()) {
    GuiMainWindow::pointer()->storeSurfaceProjection();
  }
  
  vtkIdType id_node = GuiMainWindow::pointer()->getPickedPoint();
  
  int bc_src = 42;
  int bc_dst = 18;
  int Npoints = grid->GetNumberOfPoints();
//     id_node= QInputDialog::getInt(NULL, "id_node","id_node:", id_node, 0, Npoints, 1, &ok1);
//     bc = QInputDialog::getInt(NULL, "bc","bc:", bc, 0, 100, 1, &ok2);
  QSet<int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
//     foreach (int bc, bcs) {
//         GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(grid);
//     }

  GuiMainWindow::pointer()->getSurfProj(bc_dst)->setForegroundGrid(grid);
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes  = getPartNodes();
  QVector <vtkIdType> cells;
  QSet <int> bc_src_set;
  bc_src_set.insert(bc_src);
  getSurfaceCells(bc_src_set,cells,grid);
  EG_VTKDCC(vtkIntArray, cell_code,   grid, "cell_code");
  
  foreach(vtkIdType id_cell, cells) {
    vtkIdType *pts, N_pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    for(int i=0;i<N_pts;i++) {
      vtkIdType id_node=pts[i];
      int i_nodes = _nodes[id_node];
      vec3_t x_old;
      grid->GetPoint(id_node, x_old.data());
      vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc_dst)->project(x_old, id_node);
      grid->GetPoints()->SetPoint(id_node, x_new.data());
    }
  }
  
/*  foreach(vtkIdType id_node, nodes) {
    if(cell_code->GetValue(id_node)==bc_src) {
    }
  }*/
  grid->Modified();
}
