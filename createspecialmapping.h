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
    m_grid->GetCellPoints(D, N_pts, pts);
    for(int i=0;i<N_pts;i++)
    {
      vec3_t A,B;
      m_grid->GetPoints()->GetPoint(pts[i], A.data());
      m_grid->GetPoints()->GetPoint(pts[(i+1)%N_pts], B.data());
      ret+=(B-A).abs();
    }
    return(ret);
  }
  double A_U(vtkIdType D) { // area of the circumscribed circle of the triangle
    vtkIdType N_pts, *pts;
    m_grid->GetCellPoints(D, N_pts, pts);
    vec3_t A,B,C;
    m_grid->GetPoints()->GetPoint(pts[0], A.data());
    m_grid->GetPoints()->GetPoint(pts[1], B.data());
    m_grid->GetPoints()->GetPoint(pts[2], C.data());
    double a=(C-B).abs();
    double alpha=angle((B-A),(C-A));
    double R=a/(2*sin(alpha));
    return(M_PI*R*R);
  }
  double A_D(vtkIdType D) { // triangle area
    return(cellVA(m_grid,D));
  }
  double DN(int i,vtkIdType D) { //triangle neighbours
    return(c2c[D][i]);
  }
  double nk(vtkIdType P) {
    return(n2n[P].size());
  }
  double G_k(vtkIdType node) {
//     return(CurrentMeshDensity(node,n2n,m_grid));
//     return(CurrentVertexAvgDist(node,n2n,m_grid));
    EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
    return(1.0/node_meshdensity->GetValue(node));
//     node_meshdensity->SetValue(node, D);
  };
  double DK(int i,vtkIdType D) { // triangle nodes
    vtkIdType N_pts, *pts;
    m_grid->GetCellPoints(D, N_pts, pts);
    return(pts[i]);
  }
  vtkIdType KK(int i,int j,vtkIdType K) {//i=1 or 2, j=edge nb, K=node
    if(i==1) return(K);
    else
    {
      QSet <int> S=n2n[K];
      QVector <int> V(S.size());
      qCopy(S.begin(),S.end(),V.begin());
      qSort(V.begin(),V.end());
      return(V[j]);
    }
  }
  double L_k(int j,vtkIdType K)
  {
      QSet <int> S=n2n[K];
      QVector <int> V(S.size());
      qCopy(S.begin(),S.end(),V.begin());
      qSort(V.begin(),V.end());
      vec3_t A;
      vec3_t B;
      m_grid->GetPoints()->GetPoint(K, A.data());
      m_grid->GetPoints()->GetPoint(V[j], B.data());
      return((B-A).abs());
  }
  double Q_L(vtkIdType D)
  {
    // Um(D)/sum(G_k(DK(i,D)),i,1,3)
    double denom_sum=0;
    for(int i=0;i<3;i++)
    {
      denom_sum += G_k(DK(i,D));
    }
    return(Um(D)/denom_sum);
  }
  double Q_L1(vtkIdType P)
  {
    // [2*sum(L_k(i~),i,1,nk(P))]/[sum(G_k(KK(1,i~))+G_k(KK(2,i~)),i,1,nk(P))]
    int N=nk(P);
    double num_sum=0;
    double denom_sum=0;
    for(int j=0;j<N;j++)
    {
      num_sum += 2*L_k(j,P);
      denom_sum += G_k(KK(1,j,P))+G_k(KK(2,j,P));
    }
    return(num_sum/denom_sum);
  }
  double Q_L2(vtkIdType P)
  {
    
    // min([2*L_k(i~)]/[G_k(KK(1,i~))+G_k(KK(2,i~))])
    int N=nk(P);
    QVector <double> V(N);
    double num,denom;
    for(int j=0;j<N;j++)
    {
      num = 2*L_k(j,P);
      denom = G_k(KK(1,j,P))+G_k(KK(2,j,P));
      V[j]=num/denom;
    }
    qSort(V.begin(),V.end());
    return(V[0]);
  }
  double T_min(int w)
  {
    // sum([A_U(i)]/[A_D(i)^w]*[G_k(i)^(2*(w-1))],i,1,Nd)
    int N_cells=m_grid->GetNumberOfCells();
    double T=0;
    for(int i=0;i<N_cells;i++)
    {
      T += A_U(i)/pow(A_D(i),w)*pow(G_k(i),2*(w-1));
    }
    return(T);
  }
  bool insert_fieldpoint(vtkIdType D)
  {
    double Fred1=1.0/sqrt(3);
    double Qmin=1.1;//1.189;
    double total=0;
    for(int i=0;i<3;i++)
    {
      vtkIdType cell=DN(i,D);
      if(cell!=-1) total += Q_L(cell);
    }
    cout<<"Q_L(D)>1.0/Fred1="<<Q_L(D)<<">"<<1.0/Fred1<<endl;
    cout<<"total>3*Qmin="<<total<<">"<<3*Qmin<<endl;
    return ( Q_L(D)>1.0/Fred1 && total>3*Qmin );
  }
  bool insert_edgepoint(int j,vtkIdType K)
  {
    return ( 0.5*G_k(K)<L_k(j,K) && L_k(j,K)<1*G_k(K) );
  }
  bool remove_fieldpoint(vtkIdType P)
  {
    double QL1max=0.8;
    double QL2max=0.5;
    cout<<"Q_L1(P)<QL1max="<< Q_L1(P)<< "<" << QL1max<<endl;
    cout<<"Q_L2(P)<QL2max="<< Q_L2(P)<< "<" << QL2max<<endl;
    return ( Q_L1(P)<QL1max && Q_L2(P)<QL2max );
  }
  bool remove_edgepoint(vtkIdType P)
  {
    return ( 0.5*G_k(P)<CurrentVertexAvgDist(P,n2n,m_grid) && CurrentVertexAvgDist(P,n2n,m_grid)<1*G_k(P) );
  }
};

#define VTK_SIMPLE_VERTEX 0
#define VTK_FIXED_VERTEX 1
#define VTK_FEATURE_EDGE_VERTEX 2
#define VTK_BOUNDARY_EDGE_VERTEX 3

const char* VertexType2Str(char T);
char Str2VertexType(QString S);

#endif
