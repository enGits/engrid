#ifndef GUIINFLATEMESH_H
#define GUIINFLATEMESH_H

class GuiInflateMesh;

#include "dialogoperation.h"
#include "surfacealgorithm.h"
#include "ui_guiinflatemesh.h"

class GuiInflateMesh : public DialogOperation<Ui::GuiInflateMesh, SurfaceOperation>
{

private:

  double m_EdgeLength;

  void setupSurfaceParameters();

public:

  GuiInflateMesh();

};

#endif // GUIINFLATEMESH_H
