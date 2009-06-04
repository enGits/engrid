//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
// +                                                                      +
// + enGrid is free software: you can redistribute it and/or modify       +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#include "removepoints.h"

RemovePoints::RemovePoints()
: SurfaceOperation()
{
  setQuickSave(true);
}

void RemovePoints::operate()
{
  if(remove_FP) remove_FP_all();
  if(remove_EP) remove_EP_all();
}

bool RemovePoints::remove_fieldpoint(vtkIdType P)
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

bool RemovePoints::remove_edgepoint(vtkIdType P)
{
  return ( 0.5*G_k(P)<CurrentVertexAvgDist(P) && CurrentVertexAvgDist(P)<1*G_k(P) );
}

int RemovePoints::remove_FP_counter()
{
  cout<<"===remove_FP_counter() START==="<<endl;
  QTime start = QTime::currentTime();
  
  int l_N_removed_FP = 0;
  cout<<"m_marked_cells="<<m_marked_cells<<endl;
//   cout<<"m_hitlist="<<m_hitlist<<endl;
  cout<<"m_hitlist.size()="<<m_hitlist.size()<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"l_N_removed_FP="<<l_N_removed_FP<<endl;
  
//   UpdateNodeType();
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
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
      if( !marked && remove_fieldpoint(node) && FindSnapPoint(grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        if(DebugLevel>1) cout<<"removing field point "<<node<<endl;
        l_N_removed_FP++;
        m_hitlist[node]=1;
        foreach(vtkIdType C,n2c[node]) m_marked_cells[C]=true;
        N_newcells-=2;
        N_newpoints-=1;
      }
    }
  }
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===remove_FP_counter() END==="<<endl;
  return(l_N_removed_FP);
}

int RemovePoints::remove_EP_counter()
{
  cout<<"===remove_EP_counter() START==="<<endl;
  QTime start = QTime::currentTime();
  
  int l_N_removed_EP = 0;
//   UpdateNodeType();
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
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
      if( !marked && remove_edgepoint(node) && FindSnapPoint(grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        if(DebugLevel>0) cout<<"removing edge point "<<node<<endl;
        l_N_removed_EP++;
        m_hitlist[node]=2;
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
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===remove_EP_counter() END==="<<endl;
  return(l_N_removed_EP);
}

//count all to remove, then remove them all at once
int RemovePoints::remove_FP_all()
{
  cout<<"===remove_FP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  m_hitlist.clear();
  m_offset.clear();
  m_hitlist.resize(N_points);
  m_offset.resize(N_points);
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  cout<<"================="<<endl;
  cout<<"m_hitlist.size()="<<m_hitlist.size()<<endl;
  cout<<"================="<<endl;
  
  QSet <vtkIdType> DeadNodes;
  for(vtkIdType i=0;i<m_hitlist.size();i++)
  {
    if(m_hitlist[i]==1) DeadNodes.insert(i);
  }
  int N_newpoints=0;
  int N_newcells=0;
  DeleteSetOfPoints(grid, DeadNodes, N_newpoints, N_newcells);
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  
  int kills=-N_newpoints;
  int contracts=DeadNodes.size();
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  //TEMPORARY!!!!!!!!!!!!
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===remove_FP_all END==="<<endl;
  return(0);
}

//count all to remove, then remove them all at once
int RemovePoints::remove_EP_all()
{
  cout<<"===remove_EP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  m_hitlist.clear();
  m_offset.clear();
  m_hitlist.resize(N_points);
  m_offset.resize(N_points);
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  cout<<"================="<<endl;
  cout<<"m_hitlist.size()="<<m_hitlist.size()<<endl;
  cout<<"================="<<endl;
  
  QSet <vtkIdType> DeadNodes;
  for(vtkIdType i=0;i<m_hitlist.size();i++)
  {
    if(m_hitlist[i]==1) DeadNodes.insert(i);
  }
  int N_newpoints=0;
  int N_newcells=0;
  DeleteSetOfPoints(grid, DeadNodes, N_newpoints, N_newcells);
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  
  int kills=-N_newpoints;
  int contracts=DeadNodes.size();
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===remove_EP_all END==="<<endl;
  return(0);
}
