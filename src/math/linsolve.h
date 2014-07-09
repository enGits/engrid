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
#ifndef linsolve_H
#define linsolve_H

struct LinSolveError {
  double det;
  LinSolveError(double d) { det = d; };
};

// Rainers full matrix solver
template <class M, class V>
void linsolve(const M &Ain, const V &rsv, V &b)
{
  // Solve linear system Ax = rsv
  // copy vec to protect it. b will be the return value
  b = rsv;
  // copy yourself to protect matrix entries
  M A = Ain;
  
  int n = A.size();
  int k,i,j,p[n];
  double q,s,max,h,det;
  double ele_max = 0;
  
  // Find maximum element to get a relative value
  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[i].size(); ++j) {
      ele_max = std::max(ele_max,A[i][j]);
    };
  };
  
  // Get in matrix reduction
  det = 1;
  for (k = 0; k < n-1; k++) {
    max = 0.0;
    p[k] = 0;
    for(i = k; i < n; i++) {
      s=0.0;
      for(j = k; j < n; j++) {
        s = s + fabs(A[i][j]);
      }
      q = fabs(A[i][k])/s;
      if(q > max) {
        max=q;
        p[k]=i;
      }
    }
    if(!(p[k] == k)) {
      det = -det;
      for(j = 0; j < n; j++) {
        h = A[k][j];
        A[k][j] = A[p[k]][j];
        A[p[k]][j] = h;
      }
    }
    det = det*A[k][k];
    for(i = k+1; i < n; i++) {
      A[i][k] = A[i][k]/A[k][k];
      for(j = k+1; j < n; j++) {
        A[i][j] = A[i][j] - A[i][k]*A[k][j];
      };
    }
  }
  det = det*A[n-1][n-1];
  
  // Proceed with rest of system reduction
  for(k = 0; k < n-1; k++) {
    if(!(p[k]==k)) {
      h=b[k];
      b[k]=b[p[k]];
      b[p[k]]=h;
    }
  };
  for(i = 0; i < n; i++) {
    for(j = 0; j < i; j++) {
      b[i] = b[i] - A[i][j]*b[j];
    };
  }
  for(i = n-1;i >= 0; i--) {
    s = b[i];
    for(k = i+1; k < n; k++)
      s = s - A[i][k]*b[k];
    b[i] = s/A[i][i];
  };
  
  // Check Determinant and throw error, if needed
  if (fabs(det) < 1e-20) {
    throw LinSolveError(det);
  };
};

#endif
