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
#include "brlcadinterface.h"

class CreateBrlCadTesselation : public CreateCadTesselation, public BrlCadInterface
{

public:

  CreateBrlCadTesselation(QString file_name, QString object_name);

  virtual bool shootRay(vec3_t x, vec3_t v, vec3_t &x_in, vec3_t &x_out, vec3_t &n_in, vec3_t &n_out);

};

#endif // CREATEBRLCADTESSELATION_H
