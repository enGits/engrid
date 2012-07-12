//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

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
  struct application  m_Ap;
  char                m_IdBuf[132];

  vec3_t m_LastNormal;

protected: // methods

  static int hit(struct application *ap, struct partition *PartHeadp, struct seg *segs);
  static int miss(register struct application *ap);

  vec3_t project(vec3_t x, vtkIdType id_node);

public:

  BrlCadProjection(QString file_name, QString object_name);
  ~BrlCadProjection();

  virtual vec3_t projectRestricted(vec3_t x, vtkIdType id_node = -1, bool correct_curvature = false);
  virtual vec3_t projectFree(vec3_t x, vtkIdType id_node = -1, bool correct_curvature = false);
  virtual double getRadius(vtkIdType id_node) { return 1e10; }
  virtual vec3_t lastProjNormal() { return m_LastNormal; }

};

#endif // BRLCADPROJECTION_H
