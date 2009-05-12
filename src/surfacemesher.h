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

class SurfaceMesher : public Operation {
  public:
    SurfaceMesher();
    void operate();
  
    int N_SmoothIterations;
    int maxiter_density;

    double Convergence_meshdensity;
  
    bool insert_FP;
    bool insert_EP;
    bool remove_FP;
    bool remove_EP;
  
    bool DoSwap;
    bool DoLaplaceSmoothing;
    
    int N_inserted_FP;
    int N_inserted_EP;
    int N_removed_FP;
    int N_removed_EP;
    
    int N_points;
    int N_cells;
    int N_newpoints;
    int N_newcells;
    int m_total_N_newpoints;
    int m_total_N_newcells;
  
    QSet<int> m_bcs;
    QVector<vtkIdType> m_AllCells;
    QVector<vtkIdType> m_SelectedCells;
    vtkUnstructuredGrid* m_grid;
    vtkIdType m_newNodeId;
    
    QMap< pair<vtkIdType,vtkIdType>, vtkIdType> edge_map;
    QVector <stencil_t> StencilVector;
    QVector <vtkIdType> m_SelectedNodes;
    QVector <vtkIdType> m_AllNodes;
  
    QVector <int> hitlist;//Elements to be terminated (0=keep alive, 1=field agent to eliminate, 2=border agent to eliminate)
    QVector <int> offset;//offset caused by terminated elements
  
    QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
    
    QMap <vtkIdType,bool> marked_cells;
    QMap <vtkIdType,bool> marked_nodes;
  
    void SetInput(QSet<int> a_bcs,vtkUnstructuredGrid* a_grid)
    {
      m_bcs=a_bcs;
      m_grid=a_grid;
    };
  
    void SetVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector){VMDvector=a_VMDvector;};
    void SetConvergence_meshdensity(double C){Convergence_meshdensity=C;};
    void Set_insert_FP(bool B){insert_FP=B;};
    void Set_insert_EP(bool B){insert_EP=B;};
    void Set_remove_FP(bool B){remove_FP=B;};
    void Set_remove_EP(bool B){remove_EP=B;};
  
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
      EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
      return(1.0/node_meshdensity->GetValue(node));
    };
    double DK(int i,vtkIdType D) { // triangle nodes
      vtkIdType N_pts, *pts;
      m_grid->GetCellPoints(D, N_pts, pts);
      return(pts[i]);
    }
    vtkIdType KK(int i,vtkIdType j,vtkIdType K) {//i=1 or 2, j=node2, K=node1
      if(i==1) return(K);
      else return(j);
    }
    double L_k(vtkIdType j,vtkIdType K)// node1 K, node2 j
    {
        vec3_t A;
        vec3_t B;
        m_grid->GetPoints()->GetPoint(K, A.data());
        m_grid->GetPoints()->GetPoint(j, B.data());
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
      /*if(DebugLevel>0) cout<<"D="<<D<<" Um(D)="<<Um(D)<<" denom_sum="<<denom_sum<<endl;*/
      return(Um(D)/denom_sum);
    }
    double Q_L1(vtkIdType P)
    {
      // [2*sum(L_k(i~),i,1,nk(P))]/[sum(G_k(KK(1,i~))+G_k(KK(2,i~)),i,1,nk(P))]
      double num_sum=0;
      double denom_sum=0;
      foreach(vtkIdType j,n2n[P])
      {
        num_sum += 2*L_k(j,P);
        denom_sum += G_k(KK(1,j,P))+G_k(KK(2,j,P));
      }
      return(num_sum/denom_sum);
    }
    double Q_L2(vtkIdType P)
    {
      
      // min([2*L_k(i~)]/[G_k(KK(1,i~))+G_k(KK(2,i~))])
      QVector <double> V;
      double num,denom;
      foreach(vtkIdType j,n2n[P])
      {
        num = 2*L_k(j,P);
        denom = G_k(KK(1,j,P))+G_k(KK(2,j,P));
        V.push_back(num/denom);
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
/*      cout<<"Q_L(D)>1.0/Fred1="<<Q_L(D)<<">"<<1.0/Fred1<<endl;
      cout<<"total>3*Qmin="<<total<<">"<<3*Qmin<<endl;*/
      return ( Q_L(D)>1.0/Fred1 && total>3*Qmin );
    }
    bool insert_edgepoint(vtkIdType j,vtkIdType K)// node1 K, node2 j
    {
/*      cout<<"j="<<j<<endl;
      cout<<"K="<<K<<endl;
      cout<<"0.5*G_k(K)="<<0.5*G_k(K)<<endl;
      cout<<"L_k(j,K)="<<L_k(j,K)<<endl;
      cout<<"1*G_k(K)="<<1*G_k(K)<<endl;
      cout<<"return ( 0.5*G_k(K)<L_k(j,K) && L_k(j,K)<1*G_k(K) );"<<endl;
      return ( 0.5*G_k(K)<L_k(j,K) && L_k(j,K)<1*G_k(K) );*/
      
      bool result=L_k(j,K)>0.5*(G_k(j)+G_k(K));
      if(DebugLevel>0 && result){
        cout<<"j="<<j<<endl;
        cout<<"K="<<K<<endl;
        cout<<"G_k(j)="<<G_k(j)<<endl;
        cout<<"G_k(K)="<<G_k(K)<<endl;
        cout<<"0.5*(G_k(j)+G_k(K))="<<0.5*(G_k(j)+G_k(K))<<endl;
        cout<<"L_k(j,K)="<<L_k(j,K)<<endl;
      }
      return ( result );
    }
    bool remove_fieldpoint(vtkIdType P)
    {
      double QL1max=0.8;
      double QL2max=0.5;
      bool result = Q_L1(P)<QL1max && Q_L2(P)<QL2max;
      if(DebugLevel>0 && result)
      {
        cout<<"Q_L1(P)<QL1max="<< Q_L1(P)<< "<" << QL1max<<endl;
        cout<<"Q_L2(P)<QL2max="<< Q_L2(P)<< "<" << QL2max<<endl;
      }
      return ( result );
    }
    bool remove_edgepoint(vtkIdType P)
    {
      return ( 0.5*G_k(P)<CurrentVertexAvgDist(P,n2n,m_grid) && CurrentVertexAvgDist(P,n2n,m_grid)<1*G_k(P) );
    }
    int insert_FP_counter();
    int insert_EP_counter();
    int remove_FP_counter();
    int remove_EP_counter();
    
    int insert_FP_actor(vtkUnstructuredGrid* grid_tmp);
    int insert_EP_actor(vtkUnstructuredGrid* grid_tmp);

    int insert_FP_all();
    int insert_EP_all();
  
    int UpdateDesiredMeshDensity();
  
    int remove_EP_all_3();
    int remove_FP_all_3();
  
    int SwapFunction();
    int SmoothFunction();
  
};
//end of SurfaceMesher class

#endif
