#ifndef VTKEGGRIDSMOOTHPOLYDATAFILTER_H
#define VTKEGGRIDSMOOTHPOLYDATAFILTER_H

#include "egvtkobject.h"
#include <vtkSmoothPolyDataFilter.h>

class vtkEgGridSmoothPolyDataFilter : public vtkSmoothPolyDataFilter, public EgVtkObject
{
public:
    vtkEgGridSmoothPolyDataFilter();

    ~vtkEgGridSmoothPolyDataFilter();
/*public: // static methods
  
    static vtkEgGridSmoothPolyDataFilter* New();*/
  
};

#endif
