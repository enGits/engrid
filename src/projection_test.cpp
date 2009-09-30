#include "projection_test.h"
#include "surfaceprojection.h"

#include <QInputDialog>

Projection_test::Projection_test() : SurfaceOperation()
{
      EG_TYPENAME;
}

void Projection_test::operate()
{
    bool ok1,ok2;
    int id_node = 0;
    int bc = 0;
//    int Npoints = grid->GetNumberOfPoints();
    int Npoints = 100;
    id_node= QInputDialog::getInt(NULL, "id_node","id_node:", id_node, 0, Npoints, 1, &ok1);
    bc = QInputDialog::getInt(NULL, "bc","bc:", bc, 0, 100, 1, &ok2);
    if (ok1 && ok2)
    {
//        QSet<int> bcs;
//        GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
//        foreach (int bc, bcs) {
//            GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(grid);
//        }
//
//        GuiMainWindow::pointer()->getSurfProj(bc)->setForegroundGrid(grid);
//        g2l_t _nodes = getPartLocalNodes();
//        int i_nodes = _nodes[id_node];
//        vec3_t x_old;
//        vec3_t x_new = GuiMainWindow::pointer()->getSurfProj(bc)->project(x_old, id_node);
//        grid->GetPoints()->SetPoint(id_node, x_new.data());
    }
}
