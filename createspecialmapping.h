#ifndef CREATESPECIALMAPPING_H
#define CREATESPECIALMAPPING_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>

class CreateSpecialMapping{
public:
    CreateSpecialMapping();
    int Process(vtkPolyData* input);
    ~CreateSpecialMapping();
  
  double Convergence;
  int NumberOfIterations;
  double RelaxationFactor;
  int FeatureEdgeSmoothing;
  double FeatureAngle;
  double EdgeAngle;
  int BoundarySmoothing;
  int GenerateErrorScalars;
  int GenerateErrorVectors;
  
};

#endif
