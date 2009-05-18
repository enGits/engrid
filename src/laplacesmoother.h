#ifndef LAPLACESMOOTHER_H
#define LAPLACESMOOTHER_H

#include "operation.h"

class LaplaceSmoother : public Operation {
public:
    LaplaceSmoother();
    ~LaplaceSmoother();
    void operate();
  bool FlippedCells(vtkIdType id_G, vec3_t P);
    
public:
  void SetInput(QSet<int> a_bcs,vtkUnstructuredGrid* a_grid)
  {
    m_bcs=a_bcs;
    m_grid=a_grid;
  };
  void SetNumberOfIterations(int N){NumberOfIterations=N;};
  
public:
  QSet<int> m_bcs;
  vtkUnstructuredGrid* m_grid;
  int NumberOfIterations;
  
};

#endif
