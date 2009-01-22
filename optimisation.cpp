//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#include "optimisation.h"

Optimisation::Optimisation()
{
  setDeltas(1e-6);
  F = new double**[3];
  for (int i = 0; i < 3; ++i) {
    F[i] = new double*[3];
    for (int j = 0; j < 3; ++j) {
      F[i][j] = new double[3];
    };
  };
};

void Optimisation::computeDerivatives(vec3_t x)
{
  F[0][0][0] = func(x[0]-Dx, x[1]-Dy, x[2]-Dz);
  F[1][0][0] = func(x[0],    x[1]-Dy, x[2]-Dz);
  F[2][0][0] = func(x[0]+Dx, x[1]-Dy, x[2]-Dz);
  F[0][1][0] = func(x[0]-Dx, x[1],    x[2]-Dz);
  F[1][1][0] = func(x[0],    x[1],    x[2]-Dz);
  F[2][1][0] = func(x[0]+Dx, x[1],    x[2]-Dz);
  F[0][2][0] = func(x[0]-Dx, x[1]+Dy, x[2]-Dz);
  F[1][2][0] = func(x[0],    x[1]+Dy, x[2]-Dz);
  F[2][2][0] = func(x[0]+Dx, x[1]+Dy, x[2]-Dz);
  F[0][0][1] = func(x[0]-Dx, x[1]-Dy, x[2]);
  F[1][0][1] = func(x[0],    x[1]-Dy, x[2]);
  F[2][0][1] = func(x[0]+Dx, x[1]-Dy, x[2]);
  F[0][1][1] = func(x[0]-Dx, x[1],    x[2]);
  F[1][1][1] = func(x[0],    x[1],    x[2]);
  F[2][1][1] = func(x[0]+Dx, x[1],    x[2]);
  F[0][2][1] = func(x[0]-Dx, x[1]+Dy, x[2]);
  F[1][2][1] = func(x[0],    x[1]+Dy, x[2]);
  F[2][2][1] = func(x[0]+Dx, x[1]+Dy, x[2]);
  F[0][0][2] = func(x[0]-Dx, x[1]-Dy, x[2]+Dz);
  F[1][0][2] = func(x[0],    x[1]-Dy, x[2]+Dz);
  F[2][0][2] = func(x[0]+Dx, x[1]-Dy, x[2]+Dz);
  F[0][1][2] = func(x[0]-Dx, x[1],    x[2]+Dz);
  F[1][1][2] = func(x[0],    x[1],    x[2]+Dz);
  F[2][1][2] = func(x[0]+Dx, x[1],    x[2]+Dz);
  F[0][2][2] = func(x[0]-Dx, x[1]+Dy, x[2]+Dz);
  F[1][2][2] = func(x[0],    x[1]+Dy, x[2]+Dz);
  F[2][2][2] = func(x[0]+Dx, x[1]+Dy, x[2]+Dz);
  
  grad_f[0] = (F[2][1][1]-F[0][1][1])/(2*Dx);
  grad_f[1] = (F[1][2][1]-F[1][0][1])/(2*Dy);
  grad_f[2] = (F[1][1][2]-F[1][1][0])/(2*Dz);
  
  J[0][0] = (F[0][1][1]-2*F[1][1][1]+F[2][1][1])/(Dx*Dx);
  J[1][1] = (F[1][0][1]-2*F[1][1][1]+F[1][2][1])/(Dy*Dy);
  J[2][2] = (F[1][1][0]-2*F[1][1][1]+F[1][1][2])/(Dz*Dz);
  
  J[0][1] = ((F[2][2][1]-F[0][2][1]) - (F[2][0][1]-F[0][0][1]))/(4*Dx*Dy);
  J[0][2] = ((F[2][1][2]-F[0][1][2]) - (F[2][1][0]-F[0][1][0]))/(4*Dx*Dz);
  J[1][2] = ((F[1][2][2]-F[1][2][2]) - (F[1][0][0]-F[1][0][0]))/(4*Dy*Dz);
  J[1][0] = J[0][1];
  J[2][0] = J[0][2];
  J[2][1] = J[1][2];
};

vec3_t Optimisation::optimise(vec3_t x) 
{ 
  computeDerivatives(x);
  mat3_t M = J.inverse();
  return (-1)*(M*grad_f);
};
