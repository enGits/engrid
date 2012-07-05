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
#ifndef CREATEBRLCADTESSELATION_H
#define CREATEBRLCADTESSELATION_H

#include "createcadtesselation.h"

#include "brlcad/vmath.h"
#include "brlcad/raytrace.h"

class CreateBrlCadTesselation : public CreateCadTesselation
{

private: // attributes

  struct application  m_Ap;
  struct rt_i        *m_Rtip;
  char                m_IdBuf[132];

protected: // methods

  virtual void shootRay(vec3_t x, vec3_t v);

public:

  CreateBrlCadTesselation(QString file_name, QString object_name);

};

#endif // CREATEBRLCADTESSELATION_H
