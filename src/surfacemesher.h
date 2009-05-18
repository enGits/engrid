#ifndef CREATESPECIALMAPPING_H
#define CREATESPECIALMAPPING_H

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCharArray.h>

#include <QSet>
#include <QVector>
#include <QString>
#include <QTextStream>
#include <QTime>

#include "egvtkobject.h"
#include "operation.h"
#include "vertexmeshdensity.h"
#include "smoothingutilities.h"
#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

#include <iostream>
using namespace std;

class SurfaceMesher : public Operation {
  public:
    SurfaceMesher();
    void operate();
  
    int N_SmoothIterations;
  
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
  
    vtkUnstructuredGrid* m_grid;
    
    QSet<int> m_bcs;
    QVector <vtkIdType> m_SelectedNodes;
    QVector <vtkIdType> m_AllNodes;
    QVector<vtkIdType> m_SelectedCells;
    QVector<vtkIdType> m_AllCells;
  
//     QMap <vtkIdType,bool> m_marked_cells;
//     QMap <vtkIdType,bool> m_marked_nodes;
  
    void SetInput(QSet<int> a_bcs,vtkUnstructuredGrid* a_grid)
    {
      m_bcs=a_bcs;
      m_grid=a_grid;
    };
  
    //Used for UpdateDesiredMeshDensity operation
    int MaxiterDensity;//used for UpdateDesiredMeshDensity operation
    void setMaxiterDensity(int a){MaxiterDensity=a;};
    QVector <VertexMeshDensity> VMDvector;//Vertices of Mass destruction
    void SetVertexMeshDensityVector(QVector <VertexMeshDensity> a_VMDvector){VMDvector=a_VMDvector;};
  
    double Convergence_meshdensity;
    void SetConvergence_meshdensity(double C){Convergence_meshdensity=C;};
  
    void Set_insert_FP(bool B){insert_FP=B;};
    void Set_insert_EP(bool B){insert_EP=B;};
    void Set_remove_FP(bool B){remove_FP=B;};
    void Set_remove_EP(bool B){remove_EP=B;};
  
    int SwapFunction();
    int SmoothFunction();
    void MeshDensityFunction();

  //debugging
public:
  QMap <vtkIdType,bool> m_marked_cells;
  QMap <vtkIdType,bool> m_marked_nodes;
  QMap< pair<vtkIdType,vtkIdType>, vtkIdType> m_edge_map;
  
  QVector <stencil_t> m_StencilVector;
  vtkIdType m_newNodeId;
  
  bool _insert_fieldpoint(vtkIdType D);
  bool _insert_edgepoint(vtkIdType j,vtkIdType K);// node1 K, node2 j
  int _insert_FP_counter();
  int _insert_EP_counter();
  int _insert_FP_actor(vtkUnstructuredGrid* grid_tmp);
  int _insert_EP_actor(vtkUnstructuredGrid* grid_tmp);
  int _insert_FP_all();
  int _insert_EP_all();
};
//end of SurfaceMesher class

#endif
