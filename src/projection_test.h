#ifndef PROJECTION_TEST_H
#define PROJECTION_TEST_H

#include "surfaceoperation.h"
#include "surfacealgorithm.h"

class Projection_test : public SurfaceAlgorithm
{
public:
    Projection_test();

protected: // methods
  void operate();
  void project_picked_point();
  void project_all_points();
  void Bezier_test();
  void checkInterpolationGrid();
  void Bezier_circle_test();
  void bezierFunctionTest();
  void bezierProjectionTest();
/*  Q_OBJECT;
  
protected: // methods
  
  virtual void before();
  virtual void operate();*/
  
};

#endif // PROJECTION_TEST_H
