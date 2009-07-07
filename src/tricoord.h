#ifndef TRICOORD_H
#define TRICOORD_H

#include "egvtkobject.h"

class TriCoord : public EgVtkObject
{

private: // attributes

  vtkUnstructuredGrid *m_SurfGrid;


public: // methods

  TriCoord();

};

#endif // TRICOORD_H
