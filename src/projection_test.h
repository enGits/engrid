#ifndef PROJECTION_TEST_H
#define PROJECTION_TEST_H

#include "surfaceoperation.h"

class Projection_test : public SurfaceOperation
{
public:
    Projection_test();

protected: // methods
  void operate();
  void project_picked_point();
  void project_all_points();
  
/*  Q_OBJECT;
  
protected: // methods
  
  virtual void before();
  virtual void operate();*/
  
};

#endif // PROJECTION_TEST_H
