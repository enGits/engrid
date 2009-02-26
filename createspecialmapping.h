#ifndef CREATESPECIALMAPPING_H
#define CREATESPECIALMAPPING_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <QSet>
#include <QVector>
#include "egvtkobject.h"
#include "operation.h"
#include "vertexmeshdensity.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

class CreateSpecialMapping : public Operation {
public:
    CreateSpecialMapping();
    int Process();
    void operate(){};
  
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
  QVector<vtkIdType> m_AllCells;
  QVector<vtkIdType> m_SelectedCells;
  vtkUnstructuredGrid* m_grid;
  
  double SV_value;
  double FV_value;
  double FEV_value;
  double BEV_value;
  
  QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
  
  void SetInput(QSet<int> a_bcs,vtkUnstructuredGrid* a_grid)
  {
    m_bcs=a_bcs;
    m_grid=a_grid;
  };

  void SetVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector){VMDvector=a_VMDvector;};
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
  
  VertexMeshDensity getVMD(vtkIdType node, char VertexType);

//utilities
public:
  double Um(vtkIdType D) {
    double ret=0;
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(D, N_pts, pts);
    for(int i=0;i<N_pts;i++)
    {
      vec3_t A,B;
      grid->GetPoints()->GetPoint(pts[i], A.data());
      grid->GetPoints()->GetPoint(pts[(i+1)%N_pts], B.data());
      ret+=(B-A).abs();
    }
    return(ret);
  }
  double A_U(vtkIdType D) {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(D, N_pts, pts);
    vec3_t A,B,C;
    grid->GetPoints()->GetPoint(pts[0], A.data());
    grid->GetPoints()->GetPoint(pts[1], B.data());
    grid->GetPoints()->GetPoint(pts[2], C.data());
    double a=(C-B).abs();
    double alpha=angle((B-A),(C-A));
    double R=a/(2*sin(alpha));
    return(M_PI*R*R);
  }
  double A_D(vtkIdType D) {
    return(cellVA(m_grid,D));
  }
  double DN(int i,vtkIdType D) {
    return(c2c[D][i]);
  }
  double nk(vtkIdType P) {
    return(n2n[P].size());
  }
  double G_k(vtkIdType node) {
    return(CurrentMeshDensity(node,n2n,m_grid));
  };
  double DK(int i,vtkIdType D) {
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(D, N_pts, pts);
    return(pts[i]);
  }
  double KK(int i,vtkIdType K);
  double L_k(int i);
  double Q_L(vtkIdType D);
  double Q_L1(vtkIdType D);
  double Q_L2(vtkIdType D);
  double T_min();
};

#define VTK_SIMPLE_VERTEX 0
#define VTK_FIXED_VERTEX 1
#define VTK_FEATURE_EDGE_VERTEX 2
#define VTK_BOUNDARY_EDGE_VERTEX 3

const char* VertexType2Str(char T);
char Str2VertexType(QString S);

#endif
