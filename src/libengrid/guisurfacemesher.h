#ifndef GUISURFACEMESHER_H
#define GUISURFACEMESHER_H

class GuiSurfaceMesher;

#include "dialogoperation.h"
#include "surfacemesher.h"

#include "ui_guisurfacemesher.h"

class GuiSurfaceMesher : public DialogOperation<Ui::GuiSurfaceMesher, SurfaceMesher>
{

protected: // methods

  virtual void before();
  virtual void operate();

};

#endif // GUISURFACEMESHER_H
