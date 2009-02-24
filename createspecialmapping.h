#ifndef CREATESPECIALMAPPING_H
#define CREATESPECIALMAPPING_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <QSet>
#include <QVector>
#include "egvtkobject.h"
#include "operation.h"
#include "vertexmeshdensity.h"

class CreateSpecialMapping : public Operation {
public:
    CreateSpecialMapping();
    int Process();
    void operate(){};
    ~CreateSpecialMapping();
  
  vtkPolyData* input;
  double Convergence;
  int NumberOfIterations;
  double RelaxationFactor;
  int FeatureEdgeSmoothing;
  double FeatureAngle;
  double EdgeAngle;
  int BoundarySmoothing;
  int GenerateErrorScalars;
  int GenerateErrorVectors;

  QSet<int> m_bcs;
  QVector<vtkIdType> m_cells;
  vtkUnstructuredGrid* m_grid;
  
  double SV_value;
  double FV_value;
  double FEV_value;
  double BEV_value;
  
  VertexMeshDensity VMD;//Vertices of Mass destruction
  
  void SetInput(QSet<int> a_bcs,vtkUnstructuredGrid* a_grid)
  {
    m_bcs=a_bcs;
    m_grid=a_grid;
  };

  void SetVertexMeshDensity(VertexMeshDensity a_VMD){VMD=a_VMD;};
  void SetConvergence(double C){Convergence=C;};
  void SetNumberOfIterations(int N){NumberOfIterations=N;};
  void SetRelaxationFactor(double RF){RelaxationFactor=RF;};
  void SetFeatureEdgeSmoothing(int FES){FeatureEdgeSmoothing=FES;};
  void SetFeatureAngle(double FA){FeatureAngle=FA;};
  void SetEdgeAngle(double EA){EdgeAngle=EA;};
  void SetBoundarySmoothing(int BS){BoundarySmoothing=BS;};
  void SetGenerateErrorScalars(int GES){GenerateErrorScalars=GES;};
  void SetGenerateErrorVectors(int GEV){GenerateErrorVectors=GEV;};
  
  void Set_SV_value(double V){SV_value=V;};
  void Set_FV_value(double V){FV_value=V;};
  void Set_FEV_value(double V){FEV_value=V;};
  void Set_BEV_value(double V){BEV_value=V;};
  
};

#define VTK_SIMPLE_VERTEX 0
#define VTK_FIXED_VERTEX 1
#define VTK_FEATURE_EDGE_VERTEX 2
#define VTK_BOUNDARY_EDGE_VERTEX 3

const char* VertexType2Str(char T);
char Str2VertexType(QString S);

#endif
