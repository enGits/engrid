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
  DebugLevel=1;
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
    EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
    
    m_SelectedNodes.clear();
    getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
    getNodesFromCells(m_AllCells, nodes, m_grid);
    setGrid(m_grid);
    setCells(m_AllCells);
    
    cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
    
    UpdateNodeType();
    EG_VTKDCN(vtkCharArray, node_type, m_grid, "node_type");
    
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
    
    ofstream file1,file2,file3,file4;
    
    //Method 3
    if(insert_FP) {
      insert_FP_all();
/*      file1.open ("file1.txt");
      cout_grid(file1,m_grid,true,true,true,true);
      file1.close();*/
    }
    
    if(insert_EP) {
      insert_EP_all();
/*      file2.open ("file2.txt");
      cout_grid(file2,m_grid,true,true,true,true);
      file2.close();*/
    }
    
    if(remove_FP) {
      remove_FP_all_2();
/*      file3.open ("file3.txt");
      cout_grid(file3,m_grid,true,true,true,true);
      file3.close();*/
    }
    
    if(remove_EP) {
      remove_EP_all_2();
/*      file4.open ("file4.txt");
      cout_grid(file4,m_grid,true,true,true,true);
      file4.close();*/
    }
    
    //Phase E : Delaunay swap
    if(DoSwap) {
      QSet<int> bcs_complement=complementary_bcs(m_bcs,m_grid,cells);
      cout<<"m_bcs="<<m_bcs<<endl;
      cout<<"bcs_complement="<<bcs_complement<<endl;
      
      SwapTriangles swap;
      swap.setGrid(m_grid);
      swap.setBoundaryCodes(bcs_complement);
      swap();
    }
    
      //Phase F : translate points to smooth grid
      //4 possibilities
      //vtk smooth 1
      //vtk smooth 2
      //laplacian smoothing with projection
      //Roland smoothing with projection
    
      //laplacian smoothing with projection
    if(DoLaplaceSmoothing) {
      LaplaceSmoother Lap;
      Lap.SetInput(m_bcs,m_grid);
      Lap.SetNumberOfIterations(N_SmoothIterations);
      Lap();
    }
    
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

bool CreateSpecialMapping::DeletePoint_2(vtkUnstructuredGrid *src, vtkIdType DeadNode)
{
  getAllSurfaceCells(m_AllCells,src);
  getSurfaceCells(m_bcs, m_SelectedCells, src);
  m_SelectedNodes.clear();
  getSurfaceNodes(m_bcs,m_SelectedNodes,src);
  getNodesFromCells(m_AllCells, nodes, src);
  setGrid(src);
  setCells(m_AllCells);
  
  //src grid info
  N_points=src->GetNumberOfPoints();
  N_cells=src->GetNumberOfCells();
  N_newpoints=-1;
  N_newcells=0;
  
  if(DeadNode<0 || DeadNode>=N_points)
  {
    cout<<"Warning: Point out of range: DeadNode="<<DeadNode<<" N_points="<<N_points<<endl;
    return(false);
  }
  
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  
  if(DebugLevel>10) {
    cout<<"BEFORE FINDSNAPPOINT"<<endl;
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
  vtkIdType SnapPoint=FindSnapPoint(src,DeadNode,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells);
  
  if(DebugLevel>0) cout<<"===>DeadNode="<<DeadNode<<" moving to SNAPPOINT="<<SnapPoint<<" DebugLevel="<<DebugLevel<<endl;
  if(SnapPoint<0) {cout<<"Sorry no possible SnapPoint found."<<endl; return(false);}
  
  //allocate
  if(DebugLevel>10) {
    cout<<"BEFORE ALLOCATION"<<endl;
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
  }
  N_points=src->GetNumberOfPoints();
  N_cells=src->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,dst);
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(N_points);
  
  //copy undead points
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < N_points; src_id_node++) {//loop through src points
    if(src_id_node!=DeadNode)//if the node isn't dead, copy it
    {
      vec3_t x;
      src->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(src, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
    else
    {
      if(DebugLevel>0) cout<<"src_id_node="<<src_id_node<<" dst_id_node="<<dst_id_node<<endl;
    }
  };
  if(DebugLevel>10) {
    cout<<"DeadCells="<<DeadCells<<endl;
    cout<<"MutatedCells="<<MutatedCells<<endl;
    cout<<"MutilatedCells="<<MutilatedCells<<endl;
  }
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType src_N_pts, *src_pts;
      vtkIdType dst_N_pts, *dst_pts;
      src->GetCellPoints(id_cell, src_N_pts, src_pts);
      
      vtkIdType type_cell = src->GetCellType(id_cell);
      if(DebugLevel>10) cout<<"-->id_cell="<<id_cell<<endl;
      if(DebugLevel>10) for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
//       src->GetCellPoints(id_cell, dst_N_pts, dst_pts);
      dst_N_pts=src_N_pts;
      dst_pts=new vtkIdType[dst_N_pts];
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        if(DebugLevel>10) cout<<"processing mutated cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==DeadNode) {
            if(DebugLevel>10) {
              cout<<"SnapPoint="<<SnapPoint<<endl;
              cout<<"OffSet[SnapPoint]="<<OffSet[SnapPoint]<<endl;
              cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
            }
            dst_pts[i]=SnapPoint-OffSet[SnapPoint];
          }
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
        if(DebugLevel>10) cout<<"--->dst_pts:"<<endl;
        if(DebugLevel>10) for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        
      }
      else if(MutilatedCells.contains(id_cell))//mutilated cell
      {
        if(DebugLevel>10) cout<<"processing mutilated cell "<<id_cell<<endl;
        
        if(type_cell==VTK_QUAD) {
          type_cell=VTK_TRIANGLE;
          dst_N_pts-=1;
        }
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        int j=0;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead
          //do nothing in case of DeadNode
        }
      }
      else//normal cell
      {
        if(DebugLevel>10) cout<<"processing normal cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
      if(DebugLevel>10) {
        cout<<"===Copying cell "<<id_cell<<" to "<<id_new_cell<<"==="<<endl;
        cout<<"src_pts:"<<endl;
        for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
        cout<<"dst_pts:"<<endl;
        for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        cout<<"OffSet="<<OffSet<<endl;
        cout<<"===Copying cell end==="<<endl;
      }
      delete dst_pts;
    }
  };
//   cout_grid(cout,dst,true,true,true,true);
  makeCopy(dst, src);
  return(true);
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
      GuiMainWindow::pointer()->QuickSave("pre-deleting_"+num1+"_"+num2+".vtu");
      
      if(DeletePoint_2(m_grid,i-kills))
      {
        kills++;
        cout<<"Kill successful"<<endl;
      }
      else
      {
        cout<<"Kill failed"<<endl;
        N_removed_EP--;
      }
      
      GuiMainWindow::pointer()->QuickSave("post-deleting_"+num1+"_"+num2+".vtu");
      
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
  cout<<"+++++++"<<endl;
  cout_grid(cout,m_grid,true,true,true,true);
  cout<<"+++++++"<<endl;
  
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
  
  remove_FP_counter();
  cout_grid(cout,m_grid);
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
      GuiMainWindow::pointer()->QuickSave("pre-deleting_"+num1+"_"+num2+".vtu");
      
      if(DeletePoint_2(m_grid,i-kills))
      {
        kills++;
        cout<<"Kill successful"<<endl;
      }
      else
      {
        cout<<"Kill failed"<<endl;
        N_removed_FP--;
      }
      
      GuiMainWindow::pointer()->QuickSave("post-deleting_"+num1+"_"+num2+".vtu");
      
    }
  }
  cout<<"Killed: "<<kills<<"/"<<contracts<<endl;
  if(kills!=contracts) {cout<<"MISSION FAILED"<<endl;EG_BUG;}
  cout<<"===remove_FP_all_2 END==="<<endl;
  return(0);
}
