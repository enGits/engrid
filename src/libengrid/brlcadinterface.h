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

#ifndef BRLCADINTERFACE_H
#define BRLCADINTERFACE_H

#if defined(WIN32) || defined(WIN64)
#  include <windows.h>
#  undef IGNORE
#endif

#include "brlcad/vmath.h"
#include "brlcad/raytrace.h"
#include "brlcad/common.h"

#include "engrid.h"
#include "utilities.h"

class BrlCadInterface
{

public: // data types

  enum HitType { Miss, HitIn, HitOut };
  enum PositionType { Inside, Outside, Surface };

private: // attributes

  struct application  m_Ap;
  struct rt_i        *m_Rtip;
  char                m_IdBuf[132];

  static vec3_t m_XIn;
  static vec3_t m_XOut;
  static vec3_t m_InNormal;
  static vec3_t m_OutNormal;
  static double m_InRadius;
  static double m_OutRadius;
  static bool   m_Hit;


private: // methods

  bool brlCadShootRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out, double &r_in, double &r_out);


protected: // methods

  static int hit(struct application *ap, struct partition *PartHeadp, struct seg *segs);
  static int miss(register struct application *ap);

  HitType shootRay(vec3_t x, vec3_t v, vec3_t &x_hit, vec3_t &n_hit, double &r);
  void setupBrlCad(QString file_name, QString object_name);
  PositionType position(vec3_t x, vec3_t n);


public:

  BrlCadInterface();

};

#endif // BRLCADINTERFACE_H
