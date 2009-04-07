#include "createspecialmapping.h"

#include <QString>
#include <QTextStream>
#include <vtkCharArray.h>

#include "smoothingutilities.h"

#include "swaptriangles.h"
#include "laplacesmoother.h"
#include "guimainwindow.h"

#include <iostream>
using namespace std;

CreateSpecialMapping::CreateSpecialMapping()
{
  DebugLevel=0;
}

int CreateSpecialMapping::Process()
{
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
    
    m_SelectedNodes.clear();
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
    
    //Method 1
//     FullEdit();
    
    //Method 2
/*    if(insert_FP) insert_FP_all();
    if(insert_EP) insert_EP_all();
    if(remove_FP) remove_FP_all();
    if(remove_EP) remove_EP_all();*/
    
    //Method 3
    if(insert_FP) {
      UpdateDesiredMeshDensity();
      insert_FP_all();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(insert_EP) {
      UpdateDesiredMeshDensity();
      insert_EP_all();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(remove_FP) {
      UpdateDesiredMeshDensity();
      remove_FP_all_2();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(remove_EP) {
      UpdateDesiredMeshDensity();
      remove_EP_all_2();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
/*    if(DoSwap) SwapFunction();
    if(DoLaplaceSmoothing) SmoothFunction();*/
    
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
    
    if(m_total_N_newpoints==0 && m_total_N_newcells==0) break;
    
  }
  
  cout<<"i_iter/NumberOfIterations="<<i_iter<<"/"<<NumberOfIterations<<endl;
  UpdateMeshDensity();
  return 1;
}
//end of process

int CreateSpecialMapping::UpdateDesiredMeshDensity()
{
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  
  m_SelectedNodes.clear();
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  setGrid(m_grid);
  setCells(m_AllCells);
  
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  UpdateNodeType();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  
  //Phase A : Calculate current mesh density
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
  }
  
    //Phase B : define desired mesh density
  cout<<"===Phase B==="<<endl;
  double diff=Convergence_meshdensity+1;
  if(DebugLevel>3) cout<<"before loop: diff="<<diff<<endl;
  bool first=true;
  int iter=0;
  int maxiter=100;
  do {
    if(DebugLevel>2) cout<<"--->diff="<<diff<<endl;
    first=true;
    foreach(vtkIdType node,m_SelectedNodes)
    {
      if(DebugLevel>2) cout<<"======>"<<endl;
      VertexMeshDensity nodeVMD = getVMD(node,node_type->GetValue(node));
      int idx=VMDvector.indexOf(nodeVMD);
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
  } while(diff>Convergence_meshdensity && !first && iter<maxiter);// if first=true, it means no new mesh density has been defined (all densities specified)
  cout<<"iter="<<iter<<endl;
  if(iter>=maxiter) cout<<"WARNING: Desired convergence factor has not been reached!"<<endl;
  return(0);
}

int CreateSpecialMapping::SwapFunction()
{
  //Phase E : Delaunay swap
  QSet<int> bcs_complement=complementary_bcs(m_bcs,m_grid,cells);
  cout<<"m_bcs="<<m_bcs<<endl;
  cout<<"bcs_complement="<<bcs_complement<<endl;
  
  SwapTriangles swap;
  swap.setGrid(m_grid);
  swap.setBoundaryCodes(bcs_complement);
  swap();
  return(0);
}

int CreateSpecialMapping::SmoothFunction()
{
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
  return(0);
}

VertexMeshDensity CreateSpecialMapping::getVMD(vtkIdType node, char VertexType)
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
  }
  VMD.BClist.resize(bc.size());
  qCopy(bc.begin(),bc.end(),VMD.BClist.begin());
  qSort(VMD.BClist.begin(),VMD.BClist.end());
  return(VMD);
}

int CreateSpecialMapping::insert_FP_counter()
{
  cout<<"===insert_FP_counter() START==="<<endl;
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
    if( !marked_cells[id_cell] && insert_fieldpoint(id_cell) )
    {
      cout<<"inserting a field point "<<id_cell<<endl;
      N_inserted_FP++;
      marked_cells[id_cell]=true;
      N_newcells+=2;
      N_newpoints+=1;
    }
  }
  cout<<"===insert_FP_counter() END==="<<endl;
  return(0);
}

int CreateSpecialMapping::insert_EP_counter()
{
  cout<<"===insert_EP_counter() START==="<<endl;
  
  //Phase C: Prepare edge_map
  cout<<"===Phase C==="<<endl;
  edge_map.clear();
  vtkIdType edgeId=1;
  foreach(vtkIdType node1,m_SelectedNodes)
  {
//       cout<<"node1="<<node1<<endl;
    foreach(vtkIdType node2,n2n[node1])
    {
      if(edge_map[OrderedPair(node1,node2)]==0) { //this edge hasn't been numbered yet
        edge_map[OrderedPair(node1,node2)]=edgeId;edgeId++;
      }
    }
  }
  cout<<"edge_map.size()="<<edge_map.size()<<endl;
  
  StencilVector.clear();
  QMapIterator< pair<vtkIdType,vtkIdType>, vtkIdType> edge_map_iter(edge_map);
      //rewind the iterator
  edge_map_iter.toFront ();
      //start loop
  while (edge_map_iter.hasNext()) {
    edge_map_iter.next();
    vtkIdType node1=edge_map_iter.key().first;
    vtkIdType node2=edge_map_iter.key().second;
    if(DebugLevel>10) cout << "--->(" << node1 << "," << node2 << ")" << ": " << edge_map_iter.value() << endl;
    QSet <int> stencil_cells_set;
    QVector <int> stencil_cells_vector;
    stencil_cells_set=n2c[node1];
    stencil_cells_set.intersect(n2c[node2]);
    if(DebugLevel>10) cout<<"stencil_cells_set="<<stencil_cells_set<<endl;
    
    stencil_cells_vector.resize(stencil_cells_set.size());
    qCopy(stencil_cells_set.begin(),stencil_cells_set.end(),stencil_cells_vector.begin());
    if(DebugLevel>10) cout<<"stencil_cells_vector="<<stencil_cells_vector<<endl;
    
    vtkIdType id_cell=stencil_cells_vector[0];
    int SideToSplit = getSide(id_cell,m_grid,node1,node2);
    if(DebugLevel>10) cout<<"SideToSplit="<<SideToSplit<<endl;
    if(DebugLevel>10) cout<<"c2c[id_cell][SideToSplit]="<<c2c[id_cell][SideToSplit]<<endl;
    if(DebugLevel>10) for(int i=0;i<3;i++) cout<<"c2c[id_cell]["<<i<<"]="<<c2c[id_cell][i]<<endl;
    stencil_t S=getStencil(id_cell,SideToSplit);
    
    bool stencil_marked=false;
    foreach(vtkIdType C,stencil_cells_vector)
    {
      if(marked_cells[C]) stencil_marked=true;
    }
    if(DebugLevel>10) cout<<"stencil_marked="<<stencil_marked<<endl;
    if(DebugLevel>10) cout<<"insert_edgepoint(node1,node2)="<<insert_edgepoint(node1,node2)<<endl;
    
    if( !stencil_marked && insert_edgepoint(node1,node2) )
    {
      if(DebugLevel>1) cout<<"inserting an edge point "<< "(" << node1 << "," << node2 << ")" << ": " << edge_map_iter.value() << endl;
      N_inserted_EP++;
      foreach(vtkIdType C,stencil_cells_vector) marked_cells[C]=true;
      StencilVector.push_back(S);
      
      if(stencil_cells_vector.size()==2)//2 cells around the edge
      {
        N_newcells+=2;
        N_newpoints+=1;
      }
      else//1 cell around the edge
      {
        N_newcells+=1;
        N_newpoints+=1;
      }
    }
    if(DebugLevel>10) cout <<"--->end of edge processing"<<endl;
  }
  cout<<"===insert_EP_counter() END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_FP_counter()
{
  cout<<"===remove_FP_counter() START==="<<endl;
  cout<<"marked_cells="<<marked_cells<<endl;
  cout<<"hitlist="<<hitlist<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_removed_FP="<<N_removed_FP<<endl;
  
  UpdateNodeType();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(node_type->GetValue(node)==VTK_SIMPLE_VERTEX)
    {
      bool marked=false;
      foreach(vtkIdType C,n2c[node])
      {
        if(marked_cells[C]) marked=true;
      }
      
      QSet <vtkIdType> DeadCells;
      QSet <vtkIdType> MutatedCells;
      QSet <vtkIdType> MutilatedCells;
      if( !marked && remove_fieldpoint(node) && FindSnapPoint(m_grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        if(DebugLevel>1) cout<<"removing field point "<<node<<endl;
        N_removed_FP++;
        hitlist[node]=1;
        foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
        N_newcells-=2;
        N_newpoints-=1;
      }
    }
  }
  cout<<"===remove_FP_counter() END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_EP_counter()
{
  cout<<"===remove_EP_counter() START==="<<endl;
  UpdateNodeType();
  EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(node_type->GetValue(node)==VTK_BOUNDARY_EDGE_VERTEX)
    {
      bool marked=false;
      foreach(vtkIdType C,n2c[node])
      {
        if(marked_cells[C]) marked=true;
      }
      QSet <vtkIdType> DeadCells;
      QSet <vtkIdType> MutatedCells;
      QSet <vtkIdType> MutilatedCells;
      if( !marked && remove_edgepoint(node) && FindSnapPoint(m_grid,node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        cout<<"removing edge point "<<node<<endl;
        N_removed_EP++;
        hitlist[node]=2;
        foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
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

int CreateSpecialMapping::insert_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_FP_actor START==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
/*    if(marked_cells[id_cell]) cout<<"--->marked_cells["<<id_cell<<"]=TRUE"<<endl;
    else cout<<"--->marked_cells["<<id_cell<<"]=FALSE"<<endl;*/
    
    if( !marked_cells[id_cell] && insert_fieldpoint(id_cell) )
    {
      cout<<"inserting a field point "<<id_cell<<endl;
      vtkIdType newBC=cell_code_tmp->GetValue(id_cell);
      cout<<"id_cell="<<id_cell<<" newBC="<<newBC<<endl;
      
      vtkIdType N_pts, *pts;
      m_grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t C(0,0,0);
      
      int N_neighbours=N_pts;
      cout<<"N_neighbours="<<N_neighbours<<endl;
      vec3_t corner[4];
      vtkIdType pts_triangle[4][3];
      for(int i=0;i<N_neighbours;i++)
      {
        m_grid->GetPoints()->GetPoint(pts[i], corner[i].data());
        C+=corner[i];
      }
      C=(1/(double)N_neighbours)*C;
      addPoint(grid_tmp,m_newNodeId,C.data());
      vtkIdType intmidpoint=m_newNodeId;
      m_newNodeId++;
      
      for(int i=0;i<N_neighbours;i++)
      {
        pts_triangle[i][0]=pts[i];
        pts_triangle[i][1]=pts[(i+1)%N_neighbours];
        pts_triangle[i][2]=intmidpoint;
        if(i==0)
        {
          grid_tmp->ReplaceCell(id_cell , 3, pts_triangle[0]);
          cell_code_tmp->SetValue(id_cell, newBC);
        }
        else
        {
          vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[i]);
          cell_code_tmp->SetValue(newCellId, newBC);
        }
      }
      
    }
  }
  cout<<"===insert_FP_actor END==="<<endl;
  return(0);
}

int CreateSpecialMapping::insert_EP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_EP_actor START==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(stencil_t S,StencilVector)
  {
    if(DebugLevel>10) cout<<"S="<<S<<endl;
    vec3_t A,B;
    grid_tmp->GetPoint(S.p[1],A.data());
    grid_tmp->GetPoint(S.p[3],B.data());
    vec3_t M=0.5*(A+B);
    addPoint(grid_tmp,m_newNodeId,M.data());
    
    vtkIdType pts_triangle[4][3];
    
    if(S.valid){//there is a neighbour cell
      if(DebugLevel>10) cout<<"marked_cells["<<S.id_cell1<<"]=true;"<<endl;
      if(DebugLevel>10) cout<<"marked_cells["<<S.id_cell2<<"]=true;"<<endl;
      marked_cells[S.id_cell1]=true;
      marked_cells[S.id_cell2]=true;
      
      for(int i=0;i<4;i++)
      {
        pts_triangle[i][0]=S.p[i];
        pts_triangle[i][1]=S.p[(i+1)%4];
        pts_triangle[i][2]=m_newNodeId;
      }
      
      int bc1=cell_code_tmp->GetValue(S.id_cell1);
      int bc2=cell_code_tmp->GetValue(S.id_cell2);
      
      grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
      cell_code_tmp->SetValue(S.id_cell1, bc1);
      
      grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
      cell_code_tmp->SetValue(S.id_cell2, bc2);
      
      vtkIdType newCellId;
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
      cell_code_tmp->SetValue(newCellId, bc2);
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
      cell_code_tmp->SetValue(newCellId, bc1);
    }
    else{//there is no neighbour cell
      if(DebugLevel>10) cout<<"marked_cells["<<S.id_cell1<<"]=true;"<<endl;
      marked_cells[S.id_cell1]=true;
      
      pts_triangle[0][0]=S.p[0];
      pts_triangle[0][1]=S.p[1];
      pts_triangle[0][2]=m_newNodeId;
      pts_triangle[3][0]=S.p[3];
      pts_triangle[3][1]=S.p[0];
      pts_triangle[3][2]=m_newNodeId;
      
      int bc1=cell_code_tmp->GetValue(S.id_cell1);
      
      grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
      cell_code_tmp->SetValue(S.id_cell1, bc1);
      
      vtkIdType newCellId;
      newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
      cell_code_tmp->SetValue(newCellId, bc1);
    }
    
    m_newNodeId++;
  }
  cout<<"===insert_EP_actor END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===remove_FP_actor START==="<<endl;
  
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(hitlist[node]==1)
    {
    
    }
    bool marked=false;
    foreach(vtkIdType C,n2c[node])
    {
      if(marked_cells[C]) marked=true;
    }
    if( !marked && remove_fieldpoint(node) )
    {
      if(DebugLevel>1) cout<<"removing field point "<<node<<endl;
      foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
      //TODO: Special copy function, leaving out nodes to remove
    }
  }
  cout<<"===remove_FP_actor END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_EP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===remove_EP_actor START==="<<endl;
  
  foreach(vtkIdType node,m_SelectedNodes)
  {
    bool marked=false;
    foreach(vtkIdType C,n2c[node])
    {
      if(marked_cells[C]) marked=true;
    }
    if( !marked && remove_edgepoint(node) )
    {
      cout<<"removing edge point "<<node<<endl;
      foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
      if(n2n[node].size()==4)//4 cells around the edge
      {
        
      }
      else//2 cells around the edge
      {
        
      }
    }
  }
  cout<<"===remove_EP_actor END==="<<endl;
  return(0);
}

int CreateSpecialMapping::insert_FP_all()
{
  cout<<"===insert_FP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  setGrid(m_grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_inserted_FP=0;
  
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  marked_cells.clear();
  marked_nodes.clear();
  
  insert_FP_counter();
  
    //unmark cells (TODO: optimize)
  marked_cells.clear();
    //init grid_tmp
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  insert_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  cout<<"===insert_FP_all END==="<<endl;
  return(0);
}

int CreateSpecialMapping::insert_EP_all()
{
  cout<<"===insert_EP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
  getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
  getNodesFromCells(m_AllCells, nodes, m_grid);
  setGrid(m_grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_inserted_EP=0;
  
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  marked_cells.clear();
  marked_nodes.clear();
  
  insert_EP_counter();
  
    //unmark cells (TODO: optimize)
  marked_cells.clear();
    //init grid_tmp
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  insert_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  
  cout<<"===insert_EP_all END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_FP_all()
{
  cout<<"===remove_FP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
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
  
  hitlist.resize(N_points);
  offset.resize(N_points);
  
  marked_cells.clear();
  marked_nodes.clear();
  
  remove_FP_counter();
  
    //unmark cells (TODO: optimize)
  marked_cells.clear();
    //init grid_tmp
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  remove_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  
  cout<<"===remove_FP_all END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_EP_all()
{
  cout<<"===remove_EP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
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
  
  hitlist.resize(N_points);
  offset.resize(N_points);
  
  marked_cells.clear();
  marked_nodes.clear();
  
  remove_EP_counter();
  
    //unmark cells (TODO: optimize)
  marked_cells.clear();
    //init grid_tmp
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  remove_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  
  cout<<"===remove_EP_all END==="<<endl;
  return(0);
}

int CreateSpecialMapping::FullEdit()
{
  cout<<"===FullEdit START==="<<endl;
  
  N_inserted_FP=0;
  N_inserted_EP=0;
  N_removed_FP=0;
  N_removed_EP=0;
  
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  hitlist.resize(N_points);
  offset.resize(N_points);
  
  marked_cells.clear();
  marked_nodes.clear();
  
  if(insert_FP) insert_FP_counter();
  if(insert_EP) insert_EP_counter();
  if(remove_FP) remove_FP_counter();
  if(remove_EP) remove_EP_counter();
  
  cout<<"================="<<endl;
  cout<<"hitlist="<<hitlist<<endl;
  cout<<"================="<<endl;
  
    //unmark cells (TODO: optimize)
  marked_cells.clear();
    //init grid_tmp
  N_points=m_grid->GetNumberOfPoints();
  N_cells=m_grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(m_grid, grid_tmp);//TODO: This will not work if the size of the grid is reduced!
    //initialize new node counter
  m_newNodeId=N_points;
  
  if(insert_FP) insert_FP_actor(grid_tmp);
  if(insert_EP) insert_EP_actor(grid_tmp);
  
  cout<<"================="<<endl;
  cout<<"hitlist="<<hitlist<<endl;
  cout<<"================="<<endl;
  if(remove_FP) remove_FP_actor(grid_tmp);
  if(remove_EP) remove_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  
  cout<<"===FullEdit END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_EP_all_2()
{
  cout<<"===remove_EP_all_2 START==="<<endl;
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
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
  
  marked_cells.clear();
  marked_nodes.clear();
  
  remove_EP_counter();
  cout<<"================="<<endl;
  cout<<"hitlist="<<hitlist<<endl;
  cout<<"================="<<endl;
  
  int kills=0;
  int contracts=0;
  for(int i=0;i<hitlist.size();i++)
  {
    if(hitlist[i]==2){
      contracts++;
      cout<<"Deleting point "<<i<<" currently known as "<<i-kills<<endl;
      
      QString num1;num1.setNum(i);
      QString num2;num2.setNum(i-kills);
//       GuiMainWindow::pointer()->QuickSave("pre-deleting_"+num1+"_"+num2+".vtu");
      
      bool DelResult=DeletePoint_2(m_grid,i-kills,N_newpoints,N_newcells);
      m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
      
      if(DelResult)
      {
        kills++;
        cout<<"Kill successful"<<endl;
      }
      else
      {
        cout<<"Kill failed"<<endl;
        N_removed_EP--;
      }
      
//       GuiMainWindow::pointer()->QuickSave("post-deleting_"+num1+"_"+num2+".vtu");
      
    }
  }
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  cout<<"===remove_EP_all_2 END==="<<endl;
  return(0);
}

int CreateSpecialMapping::remove_FP_all_2()
{
  cout<<"===remove_FP_all_2 START==="<<endl;
/*  cout<<"+++++++"<<endl;
  cout_grid(cout,m_grid,true,true,true,true);
  cout<<"+++++++"<<endl;*/
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
  m_SelectedNodes.clear();
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
  
  marked_cells.clear();
  marked_nodes.clear();
  
//   DualSave("pre-counter");
  remove_FP_counter();
//   DualSave("post-counter");
  
//   cout_grid(cout,m_grid);
  cout<<"================="<<endl;
  cout<<"hitlist="<<hitlist<<endl;
  cout<<"================="<<endl;
  
  int kills=0;
  int contracts=0;
  for(int i=0;i<hitlist.size();i++)
  {
    if(hitlist[i]==1){
      contracts++;
      cout<<"Deleting point "<<i<<" currently known as "<<i-kills<<endl;
      
      QString num1;num1.setNum(i);
      QString num2;num2.setNum(i-kills);
//       GuiMainWindow::pointer()->QuickSave("pre-deleting_"+num1+"_"+num2+".vtu");
      
      bool DelResult=DeletePoint_2(m_grid,i-kills,N_newpoints,N_newcells);
      m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
      
      if(DelResult)
      {
        kills++;
        cout<<"Kill successful"<<endl;
      }
      else
      {
        cout<<"Kill failed"<<endl;
        N_removed_FP--;
      }
      
//       GuiMainWindow::pointer()->QuickSave("post-deleting_"+num1+"_"+num2+".vtu");
      
    }
  }
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  cout<<"===remove_FP_all_2 END==="<<endl;
  return(0);
}
