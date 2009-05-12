//
// C++ Implementation: updatedesiredmeshdensity
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "updatedesiredmeshdensity.h"

#include <vtkCharArray.h>

UpdateDesiredMeshDensity::UpdateDesiredMeshDensity()
 : Operation()
{
}


UpdateDesiredMeshDensity::~UpdateDesiredMeshDensity()
{
}

void UpdateDesiredMeshDensity::operate()
{
  //Phase B : define desired mesh density
  cout<<"=== UpdateDesiredMeshDensity ==="<<endl;
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  getNodesFromCells(m_AllCells, m_AllNodes, grid);
  
  setGrid(grid);
  setCells(m_AllCells);
  
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  UpdateNodeType_all();
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");
  
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
      double L=CurrentVertexAvgDist(node,n2n,grid);
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
        double D=DesiredMeshDensity(node,n2n,grid);
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
  } while(diff>Convergence_meshdensity && !first && iter<MaxiterDensity);// if first=true, it means no new mesh density has been defined (all densities specified)
  cout<<"iter="<<iter<<endl;
  if(iter>=MaxiterDensity) cout<<"WARNING: Desired convergence factor has not been reached!"<<endl;
}

VertexMeshDensity UpdateDesiredMeshDensity::getVMD(vtkIdType node, char VertexType)
{
  VertexMeshDensity VMD;
  VMD.type=VertexType;
  VMD.density=0;
  VMD.CurrentNode=node;
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
/*  createNodeMapping(nodes, _nodes, grid);
  createNodeToCell(m_AllCells, nodes, _nodes, n2c, grid);*/
  
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
