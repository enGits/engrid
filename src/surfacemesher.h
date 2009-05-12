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
      return ( Q_L(D)>1.0/Fred1 && total>3*Qmin );
    }
    bool insert_edgepoint(vtkIdType j,vtkIdType K)// node1 K, node2 j
    {
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
  
    int UpdateDesiredMeshDensity();
  
    int insert_FP_counter();
    int insert_EP_counter();
    int insert_FP_actor(vtkUnstructuredGrid* grid_tmp);
    int insert_EP_actor(vtkUnstructuredGrid* grid_tmp);
    int insert_FP_all();
    int insert_EP_all();
  
    int remove_FP_counter();
    int remove_EP_counter();
    int remove_EP_all();
    int remove_FP_all();
  
    int SwapFunction();
    int SmoothFunction();
  
};
//end of SurfaceMesher class

#endif
