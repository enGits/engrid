#include "surfacemesher.h"

#include <QString>
#include <QTextStream>
#include <QTime>

#include <vtkCharArray.h>

#include "smoothingutilities.h"
#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"

#include <iostream>
using namespace std;

SurfaceMesher::SurfaceMesher()
{
   DebugLevel=0;
}

void SurfaceMesher::operate()
{
  QTime start = QTime::currentTime();
  
  int i_iter=0;
  for(i_iter=0;i_iter<NumberOfIterations;i_iter++)//TODO:Optimize this loop
  {
    cout<<"===ITERATION NB "<<i_iter<<"/"<<NumberOfIterations<<"==="<<endl;
    
    m_total_N_newpoints=0;
    m_total_N_newcells=0;
    
    getAllSurfaceCells(m_AllCells,m_grid);
    getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
    cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
    
    EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
    
    getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
    getNodesFromCells(m_AllCells, nodes, m_grid);
    setGrid(m_grid);
    setCells(m_AllCells);
    
    cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
    
    //Phase D: edit points
    cout<<"===Phase D==="<<endl;
    N_inserted_FP=0;
    N_inserted_EP=0;
    N_removed_FP=0;
    N_removed_EP=0;
    
    bool DEBUG=false;
    
    //Method 3
/*    if(insert_FP) {
      UpdateDesiredMeshDensity();
      insert_FP_all();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_FP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_FP-post-swap-1");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_FP-post-laplace");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_FP-post-swap-2");
    }
    if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/post-insert_FP");
    
    if(insert_EP) {
      UpdateDesiredMeshDensity();
      insert_EP_all();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_EP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_EP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/insert_EP-post-laplace");
    }
    if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/post-insert_EP");
    
    if(remove_FP) {
      UpdateDesiredMeshDensity();
      remove_FP_all();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/remove_FP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/remove_FP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/remove_FP-post-laplace");
    }
    if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/post-remove_FP");
    
    if(remove_EP) {
      UpdateDesiredMeshDensity();
      remove_EP_all();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    if(DEBUG) DualSave("/data1/home/mtaverne/Geometries/simulations/SurfaceTests/post-remove_EP");*/
    
    if(DoSwap) SwapFunction();
    if(DoLaplaceSmoothing) SmoothFunction();
    
    cout<<"===Summary==="<<endl;
    cout<<"N_inserted_FP="<<N_inserted_FP<<endl;
    cout<<"N_inserted_EP="<<N_inserted_EP<<endl;
    cout<<"N_removed_FP="<<N_removed_FP<<endl;
    cout<<"N_removed_EP="<<N_removed_EP<<endl;
    
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"m_total_N_newpoints="<<m_total_N_newpoints<<endl;
    cout<<"m_total_N_newcells="<<m_total_N_newcells<<endl;
    cout<<"============"<<endl;
    
//     if(m_total_N_newpoints==0 && m_total_N_newcells==0) break;
    if(N_inserted_FP==0 && N_inserted_EP==0 && N_removed_FP==0 && N_removed_EP==0) break;
  }
  
  cout<<"i_iter/NumberOfIterations="<<i_iter<<"/"<<NumberOfIterations<<endl;
  UpdateDesiredMeshDensity();
  UpdateMeshDensity();
  if(i_iter<NumberOfIterations) cout<<"WARNING: Exited before finishing all iterations."<<endl;
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  
}
//end of process

int SurfaceMesher::UpdateDesiredMeshDensity()
{
  //Phase B : define desired mesh density
  cout<<"=== UpdateDesiredMeshDensity ==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  getNodesFromCells(m_AllCells, m_AllNodes, m_grid);
  
  setGrid(m_grid);
  setCells(m_AllCells);
  
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  UpdateNodeType_all();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  EG_VTKDCN(vtkIntArray, node_specified_density, m_grid, "node_specified_density");
  
/*  //Phase A : Calculate current mesh density
  cout<<"===Phase A==="<<endl;
  
  foreach(vtkIdType node,m_SelectedNodes)
  {
    VertexMeshDensity nodeVMD = getVMD(node,node_type->GetValue(node));
    int idx=VMDvector.indexOf(nodeVMD);
    if(DebugLevel>3) cout<<"idx="<<idx<<endl;
    if(idx!=-1)//specified
    {
      node_meshdensity->SetValue(node, VMDvector[idx].density);
    }
    else//unspecified
    {
      double L=CurrentVertexAvgDist(node,n2n,m_grid);
      double D=1./L;
      node_meshdensity->SetValue(node, D);
    }
  }*/
  
  double diff=Convergence_meshdensity+1;
  if(DebugLevel>3) cout<<"before loop: diff="<<diff<<endl;
  bool first=true;
  int iter=0;
  do {
    if(DebugLevel>2) cout<<"--->diff="<<diff<<endl;
    first=true;
    foreach(vtkIdType node,m_AllNodes)
    {
      if(DebugLevel>2) cout<<"======>"<<endl;
      VertexMeshDensity nodeVMD = getVMD(node,node_type->GetValue(node));
      int idx=VMDvector.indexOf(nodeVMD);
      node_specified_density->SetValue(node, idx);
      if(DebugLevel>2) cout<<"------>idx="<<idx<<endl;
      if(idx!=-1)//specified
      {
        node_meshdensity->SetValue(node, VMDvector[idx].density);
      }
      else//unspecified
      {
        double D=DesiredMeshDensity(node,n2n,m_grid);
        if(first) {
          if(DebugLevel>2) {
            cout<<"------>FIRST:"<<endl;
            cout<<"------>D="<<D<<endl;
            cout<<"------>node_meshdensity->GetValue("<<node<<")="<<node_meshdensity->GetValue(node)<<endl;
            cout<<"------>D-node_meshdensity->GetValue("<<node<<")="<<D-node_meshdensity->GetValue(node)<<endl;
            cout<<"------>diff=abs(D-node_meshdensity->GetValue("<<node<<"))="<<abs(D-node_meshdensity->GetValue(node))<<endl;
          }
          diff=abs(D-node_meshdensity->GetValue(node));
          first=false;
        }
        else {
          if(DebugLevel>2) {
            cout<<"------>NOT FIRST:"<<endl;
            cout<<"------>D="<<D<<endl;
            cout<<"------>node_meshdensity->GetValue("<<node<<")="<<node_meshdensity->GetValue(node)<<endl;
            cout<<"------>D-node_meshdensity->GetValue("<<node<<")="<<D-node_meshdensity->GetValue(node)<<endl;
            cout<<"------>diff=abs(D-node_meshdensity->GetValue("<<node<<"))="<<abs(D-node_meshdensity->GetValue(node))<<endl;
            cout<<"------>diff="<<diff<<endl;
            cout<<"------>max(abs(D-node_meshdensity->GetValue("<<node<<")),diff)="<<max(abs(D-node_meshdensity->GetValue(node)),diff)<<endl;
          }
          diff=max(abs(D-node_meshdensity->GetValue(node)),diff);
        }
        node_meshdensity->SetValue(node, D);
      }
      if(DebugLevel>2) cout<<"======>"<<endl;
    }
    iter++;
  } while(diff>Convergence_meshdensity && !first && iter<maxiter_density);// if first=true, it means no new mesh density has been defined (all densities specified)
  cout<<"iter="<<iter<<endl;
  if(iter>=maxiter_density) cout<<"WARNING: Desired convergence factor has not been reached!"<<endl;
  return(0);
}

int SurfaceMesher::SwapFunction()
{
  //Phase E : Delaunay swap
  QSet<int> bcs_complement=complementary_bcs(m_bcs,m_grid,cells);
  cout<<"m_bcs="<<m_bcs<<endl;
  cout<<"bcs_complement="<<bcs_complement<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  setCells(m_AllCells);
  
  SwapTriangles swap;
  swap.setGrid(m_grid);
  swap.setBoundaryCodes(bcs_complement);
  swap();
  return(0);
}

int SurfaceMesher::SmoothFunction()
{
  cout<<"=== SmoothFunction START ==="<<endl;
  //Phase F : translate points to smooth grid
  //4 possibilities
  //vtk smooth 1
  //vtk smooth 2
  //laplacian smoothing with projection
  //Roland smoothing with projection

  //laplacian smoothing with projection
  LaplaceSmoother Lap;
  Lap.SetInput(m_bcs,m_grid);
  Lap.SetNumberOfIterations(N_SmoothIterations);
  Lap();
  cout<<"=== SmoothFunction END ==="<<endl;
  return(0);
}

VertexMeshDensity SurfaceMesher::getVMD(vtkIdType node, char VertexType)
{
  VertexMeshDensity VMD;
  VMD.type=VertexType;
  VMD.density=0;
  VMD.CurrentNode=node;
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
/*  createNodeMapping(nodes, _nodes, m_grid);
  createNodeToCell(m_AllCells, nodes, _nodes, n2c, m_grid);*/
  
  QSet <int> bc;
  foreach(vtkIdType C, n2c[node])
  {
    bc.insert(cell_code->GetValue(C));
    VMD.BCmap[cell_code->GetValue(C)]=2;
  }
  VMD.BClist.resize(bc.size());
  qCopy(bc.begin(),bc.end(),VMD.BClist.begin());
  qSort(VMD.BClist.begin(),VMD.BClist.end());
  return(VMD);
}

int SurfaceMesher::remove_FP_counter()
{
  cout<<"===remove_FP_counter() START==="<<endl;
  cout<<"m_marked_cells="<<m_marked_cells<<endl;
//   cout<<"hitlist="<<hitlist<<endl;
  cout<<"hitlist.size()="<<hitlist.size()<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_removed_FP="<<N_removed_FP<<endl;
  
  UpdateNodeType_all();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(node_type->GetValue(node)==VTK_SIMPLE_VERTEX)
    {
      bool marked=false;
      foreach(vtkIdType C,n2c[node])
      {
        if(m_marked_cells[C]) marked=true;
      }
      
      QSet <vtkIdType> DeadCells;
      QSet <vtkIdType> MutatedCells;
      QSet <vtkIdType> MutilatedCells;
      if( !marked && remove_fieldpoint(node) && FindSnapPoint(m_grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        if(DebugLevel>1) cout<<"removing field point "<<node<<endl;
        N_removed_FP++;
        hitlist[node]=1;
        foreach(vtkIdType C,n2c[node]) m_marked_cells[C]=true;
        N_newcells-=2;
        N_newpoints-=1;
      }
    }
  }
  cout<<"===remove_FP_counter() END==="<<endl;
  return(0);
}

int SurfaceMesher::remove_EP_counter()
{
  cout<<"===remove_EP_counter() START==="<<endl;
  UpdateNodeType_all();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(node_type->GetValue(node)==VTK_BOUNDARY_EDGE_VERTEX)
    {
      bool marked=false;
      foreach(vtkIdType C,n2c[node])
      {
        if(m_marked_cells[C]) marked=true;
      }
      QSet <vtkIdType> DeadCells;
      QSet <vtkIdType> MutatedCells;
      QSet <vtkIdType> MutilatedCells;
      if( !marked && remove_edgepoint(node) && FindSnapPoint(m_grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        if(DebugLevel>0) cout<<"removing edge point "<<node<<endl;
        N_removed_EP++;
        hitlist[node]=2;
        foreach(vtkIdType C,n2c[node]) m_marked_cells[C]=true;
        if(n2n[node].size()==4)//4 cells around the edge
        {
          N_newcells-=2;
          N_newpoints-=1;
        }
        else//2 cells around the edge
        {
          N_newcells-=1;
          N_newpoints-=1;
        }
      }
    }
  }
  cout<<"===remove_EP_counter() END==="<<endl;
  return(0);
}

//count all to remove, then remove them all at once
int SurfaceMesher::remove_FP_all()
{
  cout<<"===remove_FP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  setGrid(m_grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_removed_FP=0;
  
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  hitlist.clear();
  offset.clear();
  hitlist.resize(N_points);
  offset.resize(N_points);
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  remove_FP_counter();
  cout<<"================="<<endl;
  cout<<"hitlist.size()="<<hitlist.size()<<endl;
  cout<<"================="<<endl;
  
  QSet <vtkIdType> DeadNodes;
  for(vtkIdType i=0;i<hitlist.size();i++)
  {
    if(hitlist[i]==1) DeadNodes.insert(i);
  }
  int N_newpoints=0;
  int N_newcells=0;
  DeleteSetOfPoints(m_grid, DeadNodes, N_newpoints, N_newcells);
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  
  int kills=-N_newpoints;
  int contracts=DeadNodes.size();
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  cout<<"===remove_FP_all END==="<<endl;
  return(0);
}

//count all to remove, then remove them all at once
int SurfaceMesher::remove_EP_all()
{
  cout<<"===remove_EP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  setGrid(m_grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_removed_EP=0;
  
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  hitlist.clear();
  offset.clear();
  hitlist.resize(N_points);
  offset.resize(N_points);
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  remove_EP_counter();
  cout<<"================="<<endl;
  cout<<"hitlist.size()="<<hitlist.size()<<endl;
  cout<<"================="<<endl;
  
  QSet <vtkIdType> DeadNodes;
  for(vtkIdType i=0;i<hitlist.size();i++)
  {
    if(hitlist[i]==1) DeadNodes.insert(i);
  }
  int N_newpoints=0;
  int N_newcells=0;
  DeleteSetOfPoints(m_grid, DeadNodes, N_newpoints, N_newcells);
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  
  int kills=-N_newpoints;
  int contracts=DeadNodes.size();
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  cout<<"===remove_EP_all END==="<<endl;
  return(0);
}
