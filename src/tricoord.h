#ifndef TRICOORD_H
#define TRICOORD_H

#include "egvtkobject.h"

class TriCoord : public EgVtkObject
{

private: // attributes

  vtkUnstructuredGrid* m_SurfGrid;
  int m_Sector;

public: // methods

  TriCoord(vtkUnstructuredGrid* surf_grid, vtkIdType id_tri, vec3_t x);
  void setPosition(vec3_t x);

};

#endif // TRICOORD_H
