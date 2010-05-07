//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#ifndef PROJECTION_TEST_H
#define PROJECTION_TEST_H

#include "surfaceoperation.h"
#include "surfacealgorithm.h"

#include "beziertriangle.h"

class Projection_test : public SurfaceAlgorithm {
  public:
    Projection_test();

  protected: // methods
    void operate();
    void project_picked_point();
    void project_all_points();
    void project_all_points2();
    void Bezier_test();
    void checkInterpolationGrid();
    void Bezier_circle_test();
    void bezierFunctionTest();
    void bezierProjectionTest();
    void bezierQuads();
    void bezierProjectionTest2(BezierTriangle bezier_triangle, QString prefix);
    BezierTriangle specialTriangle(bool equi, int type);
    
    /*  Q_OBJECT;

    protected: // methods

      virtual void before();
      virtual void operate();*/

private: // methods
  
  bool setNewPosition(vtkIdType id_node, vec3_t x_new);
  bool moveNode(vtkIdType id_node, vec3_t &Dx);
  
};

#endif // PROJECTION_TEST_H
