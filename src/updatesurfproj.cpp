#include "updatesurfproj.h"
#include "guimainwindow.h"

UpdateSurfProj::UpdateSurfProj()
{
}

void UpdateSurfProj::operate()
{
  GuiMainWindow::pointer()->storeSurfaceProjection();
}
