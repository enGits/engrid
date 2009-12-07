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

};

#endif // PROJECTION_TEST_H
