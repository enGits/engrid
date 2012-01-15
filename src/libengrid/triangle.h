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
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <QVector>

#include "vtkIdList.h"
#include "vtkUnstructuredGrid.h"
#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "egvtkobject.h"

class Triangle : public EgVtkObject
{

protected:

  vtkIdType m_IdA;
  vtkIdType m_IdB;
  vtkIdType m_IdC;

  vec3_t m_Xa, m_Xb, m_Xc;
  vec3_t m_G1, m_G2, m_G3;
  mat3_t m_G, m_GI;
  double m_A;
  double m_SmallestLength;
  double m_SmallestHeight;
  bool   m_Valid;
  vec3_t m_NormalA, m_NormalB, m_NormalC;
  vec3_t m_RNormalA, m_RNormalB, m_RNormalC;
  QVector <bool> m_HasNeighbour; ///< True if edge i has a neighbour in the grid

public:

  Triangle();
  Triangle(vec3_t a, vec3_t b, vec3_t c);
  Triangle(vtkUnstructuredGrid* grid, vtkIdType id_a, vtkIdType id_b, vtkIdType id_c);
  Triangle(vtkUnstructuredGrid* grid, vtkIdType id_cell);
  void setupTriangle();
  void setDefaults();
  
  /**
     * Calculates the closest (NOT the projection!) point (xi,ri) of point xp on the triangle.
     * @param xp Point to "project"
     * @param xi Global 3D coordinates of the closest point on the triangle.
     * @param ri Local 3D triangle coordinates of the closest point on the triangle. (0<=ri[0]<=1 and 0<=ri[1]<=1 and ri[2]=0)
     * @param d Distance of xp to (xi,ri)
     * @return True if (xi,ri) is the result of a direct projection on the triangle, else false.
    */
  bool projectOnTriangle(vec3_t xp, vec3_t &xi, vec3_t &ri, double &d, int& side, bool restrict_to_triangle);

  vec3_t local3DToGlobal3D(vec3_t l_M);
  vec3_t global3DToLocal3D(vec3_t g_M);
  vec3_t local2DToGlobal3D(vec2_t l_M);
  vec2_t global3DToLocal2D(vec3_t g_M);

  void saveTriangle(QString filename);

  bool hasNeighbour(int i) { return m_HasNeighbour[i]; }
  void setNeighbourTrue(int i)  { m_HasNeighbour[i] = true; }
  void setNeighbourFalse(int i) { m_HasNeighbour[i] = false; }
  vtkIdType idA() { return m_IdA; }
  vtkIdType idB() { return m_IdB; }
  vtkIdType idC() { return m_IdC; }
  vec3_t g1() { return m_G1; }
  vec3_t g2() { return m_G2; }
  vec3_t g3() { return m_G3; }
  vec3_t a() { return m_Xa; }
  vec3_t b() { return m_Xb; }
  vec3_t c() { return m_Xc; }
  vec3_t nA() { return m_NormalA; }
  vec3_t nB() { return m_NormalB; }
  vec3_t nC() { return m_NormalC; }
  vec3_t rNa() { return m_RNormalA; }
  vec3_t rNb() { return m_RNormalB; }
  vec3_t rNc() { return m_RNormalC; }
  double smallestLength() { return m_SmallestLength; }
  double smallestHeight() { return m_SmallestHeight; }
  void setNormals(vec3_t na, vec3_t nb, vec3_t nc);


};

#endif
