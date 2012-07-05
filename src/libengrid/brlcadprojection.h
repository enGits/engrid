#ifndef BRLCADPROJECTION_H
#define BRLCADPROJECTION_H

class BrlCadProjection;

#include "surfaceprojection.h"

//#include "common.h"
//#include "machine.h"
#include "engrid.h"
#include "brlcad/vmath.h"
#include "brlcad/raytrace.h"

class BrlCadProjection : public SurfaceProjection
{

  struct rt_i        *m_Rtip;
  struct application  m_App;
  char                m_IdBuf[132];

  void load(QString file_name);
  void prepare(QString object_name);

public:

  BrlCadProjection(QString file_name);
  ~BrlCadProjection();

  virtual vec3_t projectRestricted(vec3_t x, vtkIdType id_node = -1);
  virtual vec3_t projectFree(vec3_t x, vtkIdType id_node = -1);

};

#endif // BRLCADPROJECTION_H
