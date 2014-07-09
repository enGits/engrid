// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
#include "elements.h"

Elements::Elements()
{
  setTet(0,0, 1,4,3,5);
  setTet(0,1, 0,3,2,1);
  setTet(0,2, 2,3,5,1);
  setTet(1,0, 1,4,3,5);
  setTet(1,1, 0,3,5,1);
  setTet(1,2, 0,5,2,1);
  setTet(2,0, 0,1,4,2);
  setTet(2,1, 0,3,5,4);
  setTet(2,2, 0,5,2,4);
  setTet(3,0, 0,1,4,2);
  setTet(3,1, 0,3,2,4);
  setTet(3,2, 3,5,2,4);
};

void Elements::setTet(int i_variant, int i_tetra, int n0, int n1, int n2, int n3)
{
  pri_tet[i_variant][i_tetra][0] = n0;
  pri_tet[i_variant][i_tetra][1] = n1;
  pri_tet[i_variant][i_tetra][2] = n2;
  pri_tet[i_variant][i_tetra][3] = n3;
};


