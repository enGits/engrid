#include "createspecialmapping.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangleFilter.h"
#include <QString>
#include <QTextStream>

#include "smoothingutilities.h"

#include "swaptriangles.h"
#include "laplacesmoother.h"

#include <iostream>
using namespace std;

const char* VertexType2Str(char T)
{
  if(T==VTK_SIMPLE_VERTEX) return("VTK_SIMPLE_VERTEX");
  if(T==VTK_FIXED_VERTEX) return("VTK_FIXED_VERTEX");
  if(T==VTK_FEATURE_EDGE_VERTEX) return("VTK_FEATURE_EDGE_VERTEX");
  if(T==VTK_BOUNDARY_EDGE_VERTEX) return("VTK_BOUNDARY_EDGE_VERTEX");
  else return("Unknown vertex type");
}

char Str2VertexType(QString S)
{
  if(S=="VTK_SIMPLE_VERTEX") return((char)0);
  if(S=="VTK_FIXED_VERTEX") return((char)1);
  if(S=="VTK_FEATURE_EDGE_VERTEX") return((char)2);
  if(S=="VTK_BOUNDARY_EDGE_VERTEX") return((char)3);
  else return((char)-1);
}

// Special structure for marking vertices
typedef struct _vtkMeshVertex 
{
  char      type;
  vtkIdList *edges; // connected edges (list of connected point ids)
} vtkMeshVertex, *vtkMeshVertexPtr;

CreateSpecialMapping::CreateSpecialMapping()
{
}

int CreateSpecialMapping::Process()
{
//   vtkMeshVertexPtr Verts;
  EG_VTKDCN(vtkIntArray, node_type, m_grid, "node_type");
  
  DebugLevel=0;
  for(int i_iter=0;i_iter<NumberOfIterations;i_iter++)//TODO:Optimize this loop
  {
    cout<<"===ITERATION NB "<<i_iter<<"/"<<NumberOfIterations<<"==="<<endl;
    getAllSurfaceCells(m_AllCells,m_grid);
    getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
    cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
    
    EG_VTKDCC(vtkIntArray, cell_code, m_grid, "cell_code");
    EG_VTKDCN(vtkDoubleArray, node_meshdensity, m_grid, "node_meshdensity");
    
    m_SelectedNodes.clear();
    getSurfaceNodes(m_bcs,m_SelectedNodes,m_grid);
    getNodesFromCells(m_AllCells, nodes, m_grid);
/*    createNodeMapping(nodes, _nodes, m_grid);
    createNodeMapping(nodes, _nodes, m_grid);
    createNodeToNode(m_AllCells, nodes, _nodes, n2n, m_grid);
    createCellToCell(m_AllCells, c2c, m_grid);*/
    setGrid(m_grid);
    setCells(m_AllCells);
    
    cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
    
    UpdateNodeType();
    
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
/*    int iter=0;
    int maxiter=100;*/
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
          double L=1./D;
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
//       iter++;
    } while(diff>Convergence_meshdensity && !first /*&& iter<maxiter*/);// if first=true, it means no new mesh density has been defined (all densities specified)
//     cout<<"iter="<<iter<<endl;
//     if(iter>=maxiter) EG_BUG;
    
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
    
    FullEdit();
    
/*    if(insert_FP) insert_FP_all();
    if(insert_EP) insert_EP_all();
    if(remove_FP) remove_FP_all();
    if(remove_EP) remove_EP_all();*/
    
    cout<<"N_inserted_FP="<<N_inserted_FP<<endl;
    cout<<"N_inserted_EP="<<N_inserted_EP<<endl;
    cout<<"N_removed_FP="<<N_removed_FP<<endl;
    cout<<"N_removed_EP="<<N_removed_EP<<endl;
    
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
    
    //Phase E : Delaunay swap
    QSet<int> bcs_complement=complementary_bcs(m_bcs,m_grid,cells);
    cout<<"m_bcs="<<m_bcs<<endl;
    cout<<"bcs_complement="<<bcs_complement<<endl;
    
    SwapTriangles swap;
    swap.setGrid(m_grid);
    swap.setBoundaryCodes(bcs_complement);
    swap();
    
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
    
    cout<<"===Summary==="<<endl;
    cout<<"N_inserted_FP="<<N_inserted_FP<<endl;
    cout<<"N_inserted_EP="<<N_inserted_EP<<endl;
    cout<<"N_removed_FP="<<N_removed_FP<<endl;
    cout<<"N_removed_EP="<<N_removed_EP<<endl;
    
    cout<<"N_points="<<N_points<<endl;
    cout<<"N_cells="<<N_cells<<endl;
    cout<<"N_newpoints="<<N_newpoints<<endl;
    cout<<"N_newcells="<<N_newcells<<endl;
    cout<<"============"<<endl;
    
  }
  
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
  createNodeMapping(nodes, _nodes, m_grid);
  createNodeToCell(m_AllCells, nodes, _nodes, n2c, m_grid);
  
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
  cout<<"===insert_FP_counter()==="<<endl;
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
  return(0);
}

int CreateSpecialMapping::insert_EP_counter()
{
  cout<<"===insert_EP_counter()==="<<endl;
  StencilVector.clear();
  QMapIterator< pair<vtkIdType,vtkIdType>, vtkIdType> edge_map_iter(edge_map);
      //rewind the iterator
  edge_map_iter.toFront ();
      //start loop
  while (edge_map_iter.hasNext()) {
    edge_map_iter.next();
    vtkIdType node1=edge_map_iter.key().first;
    vtkIdType node2=edge_map_iter.key().second;
    cout << "--->(" << node1 << "," << node2 << ")" << ": " << edge_map_iter.value() << endl;
    QSet <int> stencil_cells_set;
    QVector <int> stencil_cells_vector;
    stencil_cells_set=n2c[node1];
    stencil_cells_set.intersect(n2c[node2]);
    cout<<"stencil_cells_set="<<stencil_cells_set<<endl;
    
    stencil_cells_vector.resize(stencil_cells_set.size());
    qCopy(stencil_cells_set.begin(),stencil_cells_set.end(),stencil_cells_vector.begin());
    cout<<"stencil_cells_vector="<<stencil_cells_vector<<endl;
    
    vtkIdType id_cell=stencil_cells_vector[0];
    int SideToSplit = getSide(id_cell,m_grid,node1,node2);
    cout<<"SideToSplit="<<SideToSplit<<endl;
    cout<<"c2c[id_cell][SideToSplit]="<<c2c[id_cell][SideToSplit]<<endl;
    for(int i=0;i<3;i++) cout<<"c2c[id_cell]["<<i<<"]="<<c2c[id_cell][i]<<endl;
    stencil_t S=getStencil(id_cell,SideToSplit);
    
    bool stencil_marked=false;
    foreach(vtkIdType C,stencil_cells_vector)
    {
      if(marked_cells[C]) stencil_marked=true;
    }
    cout<<"stencil_marked="<<stencil_marked<<endl;
    cout<<"insert_edgepoint(node1,node2)="<<insert_edgepoint(node1,node2)<<endl;
    
    if( !stencil_marked && insert_edgepoint(node1,node2) )
    {
      cout<<"inserting an edge point "<< "(" << node1 << "," << node2 << ")" << ": " << edge_map_iter.value() << endl;
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
    cout <<"--->end of edge processing"<<endl;
  }
  return(0);
}

int CreateSpecialMapping::remove_FP_counter()
{
  cout<<"===remove_FP_counter()==="<<endl;
  
  foreach(vtkIdType node,m_SelectedNodes)
  {
    bool marked=false;
    foreach(vtkIdType C,n2c[node])
    {
      if(marked_cells[C]) marked=true;
    }
    if( !marked && remove_fieldpoint(node) )
    {
      cout<<"removing field point "<<node<<endl;
      N_removed_FP++;
      hitlist[node]=1;
      foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
      N_newcells-=2;
      N_newpoints-=1;
    }
  }
  return(0);
}

int CreateSpecialMapping::remove_EP_counter()
{
  cout<<"===remove_EP_counter()==="<<endl;
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
  return(0);
}

int CreateSpecialMapping::insert_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_FP_actor==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
    if(marked_cells[id_cell]) cout<<"--->marked_cells["<<id_cell<<"]=TRUE"<<endl;
    else cout<<"--->marked_cells["<<id_cell<<"]=FALSE"<<endl;
    
    if( !marked_cells[id_cell] && insert_fieldpoint(id_cell) )
    {
      cout<<"inserting a field point "<<id_cell<<endl;
      vtkIdType newBC=cell_code_tmp->GetValue(id_cell);
      cout<<"id_cell="<<id_cell<<" newBC="<<newBC<<endl;
      
      vtkIdType N_pts, *pts;
      m_grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t C(0,0,0);
      
      vtkIdType type_cell = m_grid->GetCellType(id_cell);
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
  return(0);
}

int CreateSpecialMapping::insert_EP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_EP_actor==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(stencil_t S,StencilVector)
  {
    cout<<"S="<<S<<endl;
    vec3_t A,B;
    grid_tmp->GetPoint(S.p[1],A.data());
    grid_tmp->GetPoint(S.p[3],B.data());
    vec3_t M=0.5*(A+B);
    addPoint(grid_tmp,m_newNodeId,M.data());
    
    vtkIdType pts_triangle[4][3];
    
    if(S.valid){//there is a neighbour cell
      cout<<"marked_cells["<<S.id_cell1<<"]=true;"<<endl;
      cout<<"marked_cells["<<S.id_cell2<<"]=true;"<<endl;
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
      cout<<"marked_cells["<<S.id_cell1<<"]=true;"<<endl;
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
  return(0);
}

int CreateSpecialMapping::remove_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===remove_FP_actor==="<<endl;
  
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
      cout<<"removing field point "<<node<<endl;
      foreach(vtkIdType C,n2c[node]) marked_cells[C]=true;
      //TODO: Special copy function, leaving out nodes to remove
    }
  }
  return(0);
}

int CreateSpecialMapping::remove_EP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===remove_EP_actor==="<<endl;
  
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
  return(0);
}

int CreateSpecialMapping::insert_FP_all()
{
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
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  insert_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  return(0);
}

int CreateSpecialMapping::insert_EP_all()
{
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
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  insert_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  return(0);
}

int CreateSpecialMapping::remove_FP_all()
{
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
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  remove_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  return(0);
}

int CreateSpecialMapping::remove_EP_all()
{
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
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  makeCopyNoAlloc(m_grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  remove_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,m_grid);
  return(0);
}

int CreateSpecialMapping::FullEdit()
{
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
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
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
  return(0);
}

int CreateSpecialMapping::UpdateMeshDensity()
{
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, m_grid, "node_meshdensity_current");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    double L=CurrentVertexAvgDist(node,n2n,m_grid);
    double D=1./L;
    node_meshdensity_current->SetValue(node, D);
  }
  return(0);
}

int CreateSpecialMapping::UpdateNodeType()
{
  cout<<"===UpdateNodeType==="<<endl;
  DebugLevel=6;
  
  getAllSurfaceCells(m_AllCells,m_grid);
  getSurfaceCells(m_bcs, m_SelectedCells, m_grid);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  EG_VTKSP(vtkPolyData, pdata);
  //   addToPolyData(m_SelectedCells, pdata, m_grid);
  addToPolyData(m_AllCells, pdata, m_grid);
  input=pdata;
  
  vtkPolyData *source = 0;
  
  vtkIdType numPts, numCells, i, numPolys, numStrips;
  int j, k;
  vtkIdType npts = 0;
  vtkIdType *pts = 0;
  vtkIdType p1, p2;
  double x[3], y[3], deltaX[3], xNew[3], conv, maxDist, dist, factor;
  double x1[3], x2[3], x3[3], l1[3], l2[3];
  double CosFeatureAngle; //Cosine of angle between adjacent polys
  double CosEdgeAngle; // Cosine of angle between adjacent edges
  double closestPt[3], dist2, *w = NULL;
  int iterationNumber, abortExecute;
  vtkIdType numSimple=0, numBEdges=0, numFixed=0, numFEdges=0;
  vtkPolyData *inMesh, *Mesh;
  vtkPoints *inPts;
  vtkTriangleFilter *toTris=NULL;
  vtkCellArray *inVerts, *inLines, *inPolys, *inStrips;
  vtkPoints *newPts;
  vtkMeshVertexPtr Verts;
  vtkCellLocator *cellLocator=NULL;
  
    // Check input
    //
  numPts=input->GetNumberOfPoints();
  numCells=input->GetNumberOfCells();
  if (numPts < 1 || numCells < 1)
  {
    cout<<"No data to smooth!"<<endl;
    return 1;
  }
  
  CosFeatureAngle = 
    cos((double) vtkMath::DegreesToRadians() * this->FeatureAngle);
  CosEdgeAngle = cos((double) vtkMath::DegreesToRadians() * this->EdgeAngle);
  
  if(DebugLevel>5) {
    cout<<"Smoothing " << numPts << " vertices, " << numCells 
      << " cells with:\n"
      << "\tConvergence= " << this->Convergence << "\n"
      << "\tIterations= " << this->NumberOfIterations << "\n"
      << "\tRelaxation Factor= " << this->RelaxationFactor << "\n"
      << "\tEdge Angle= " << this->EdgeAngle << "\n"
      << "\tBoundary Smoothing " << (this->BoundarySmoothing ? "On\n" : "Off\n")
      << "\tFeature Edge Smoothing " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n")
      << "\tError Scalars " << (this->GenerateErrorScalars ? "On\n" : "Off\n")
      << "\tError Vectors " << (this->GenerateErrorVectors ? "On\n" : "Off\n")<<endl;
  }
    // Peform topological analysis. What we're gonna do is build a connectivity
    // array of connected vertices. The outcome will be one of three
    // classifications for a vertex: VTK_SIMPLE_VERTEX, VTK_FIXED_VERTEX. or
    // VTK_EDGE_VERTEX. Simple vertices are smoothed using all connected 
    // vertices. FIXED vertices are never smoothed. Edge vertices are smoothed
    // using a subset of the attached vertices.
    //
  cout<<"===>Analyze topology==="<<endl;
  if(DebugLevel>5) cout<<"Analyzing topology..."<<endl;
  cout<<"0:numPts="<<numPts<<endl;
  Verts = new vtkMeshVertex[numPts];
  for (i=0; i<numPts; i++)
  {
    if(DebugLevel>5) cout<<"0:VTK_SIMPLE_VERTEX"<<endl;
    Verts[i].type = VTK_SIMPLE_VERTEX; //can smooth
    Verts[i].edges = NULL;
  }
  
  inPts = input->GetPoints();
  conv = this->Convergence * input->GetLength();
  
    // check vertices first. Vertices are never smoothed_--------------
  for (inVerts=input->GetVerts(), inVerts->InitTraversal(); 
       inVerts->GetNextCell(npts,pts); )
  {
    for (j=0; j<npts; j++)
    {
      if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"->vertices:VTK_FIXED_VERTEX"<<endl;
      Verts[pts[j]].type = VTK_FIXED_VERTEX;
    }
  }
  
  cout<<"==check lines=="<<endl;
    // now check lines. Only manifold lines can be smoothed------------
  for (inLines=input->GetLines(), inLines->InitTraversal(); 
       inLines->GetNextCell(npts,pts); )
  {
    for (j=0; j<npts; j++)
    {
      if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"->lines"<<endl;
      if ( Verts[pts[j]].type == VTK_SIMPLE_VERTEX )
      {
        if ( j == (npts-1) ) //end-of-line marked FIXED
        {
          if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"2:VTK_FIXED_VERTEX"<<endl;
          Verts[pts[j]].type = VTK_FIXED_VERTEX;
        }
        else if ( j == 0 ) //beginning-of-line marked FIXED
        {
          if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"3:VTK_FIXED_VERTEX"<<endl;
          Verts[pts[0]].type = VTK_FIXED_VERTEX;
          inPts->GetPoint(pts[0],x2);
          inPts->GetPoint(pts[1],x3);
        }
        else //is edge vertex (unless already edge vertex!)
        {
          if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"4:VTK_FEATURE_EDGE_VERTEX"<<endl;
          Verts[pts[j]].type = VTK_FEATURE_EDGE_VERTEX;
          Verts[pts[j]].edges = vtkIdList::New();
          Verts[pts[j]].edges->SetNumberOfIds(2);
          Verts[pts[j]].edges->SetId(0,pts[j-1]);
          Verts[pts[j]].edges->SetId(1,pts[j+1]);
        }
      } //if simple vertex
      
      else if ( Verts[pts[j]].type == VTK_FEATURE_EDGE_VERTEX )
      { //multiply connected, becomes fixed!
        if(DebugLevel>5) cout<<"pts[j]="<<pts[j]<<"5:VTK_FIXED_VERTEX"<<endl;
        Verts[pts[j]].type = VTK_FIXED_VERTEX;
        Verts[pts[j]].edges->Delete();
        Verts[pts[j]].edges = NULL;
      }
      
    } //for all points in this line
  } //for all lines
  
  cout<<"==polygons and triangle strips=="<<endl;
    // now polygons and triangle strips-------------------------------
  inPolys=input->GetPolys();
  numPolys = inPolys->GetNumberOfCells();
  inStrips=input->GetStrips();
  numStrips = inStrips->GetNumberOfCells();
  
  if(DebugLevel>5) cout<<"numPolys="<<numPolys<<endl;
  if(DebugLevel>5) cout<<"numStrips="<<numStrips<<endl;
  
  
  if ( numPolys > 0 || numStrips > 0 )
  { //build cell structure
    vtkCellArray *polys;
    vtkIdType cellId;
    int numNei, nei, edge;
    vtkIdType numNeiPts;
    vtkIdType *neiPts;
    double normal[3], neiNormal[3];
    vtkIdList *neighbors;
    
    neighbors = vtkIdList::New();
    neighbors->Allocate(VTK_CELL_SIZE);
    
    inMesh = vtkPolyData::New();
    inMesh->SetPoints(inPts);
    inMesh->SetPolys(inPolys);
    Mesh = inMesh;
    
    if ( (numStrips = inStrips->GetNumberOfCells()) > 0 )
    { // convert data to triangles
      inMesh->SetStrips(inStrips);
      toTris = vtkTriangleFilter::New();
      toTris->SetInput(inMesh);
      toTris->Update();
      Mesh = toTris->GetOutput();
    }
    
    Mesh->BuildLinks(); //to do neighborhood searching
    polys = Mesh->GetPolys();
    
    for (cellId=0, polys->InitTraversal(); polys->GetNextCell(npts,pts); 
         cellId++)
    {
      if(DebugLevel>5) cout<<"->cellId="<<cellId<<endl;
      for (i=0; i < npts; i++) 
      {
        if(DebugLevel>5) cout<<"-->i="<<i<<endl;
        p1 = pts[i];
        p2 = pts[(i+1)%npts];
        
        if ( Verts[p1].edges == NULL )
        {
          Verts[p1].edges = vtkIdList::New();
          Verts[p1].edges->Allocate(16,6);
        }
        if ( Verts[p2].edges == NULL )
        {
          Verts[p2].edges = vtkIdList::New();
          Verts[p2].edges->Allocate(16,6);
        }
        
        Mesh->GetCellEdgeNeighbors(cellId,p1,p2,neighbors);
        numNei = neighbors->GetNumberOfIds();
        if(DebugLevel>5) cout<<"-->numNei="<<numNei<<endl;
        
        edge = VTK_SIMPLE_VERTEX;
        if ( numNei == 0 )
        {
          edge = VTK_BOUNDARY_EDGE_VERTEX;
        }
        
        else if ( numNei >= 2 )
        {
            // check to make sure that this edge hasn't been marked already
          for (j=0; j < numNei; j++)
          {
            if ( neighbors->GetId(j) < cellId )
            {
              break;
            }
          }
          if ( j >= numNei )
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        
        else if ( numNei == 1 && (nei=neighbors->GetId(0)) > cellId ) 
        {
          vtkPolygon::ComputeNormal(inPts,npts,pts,normal);
          Mesh->GetCellPoints(nei,numNeiPts,neiPts);
          vtkPolygon::ComputeNormal(inPts,numNeiPts,neiPts,neiNormal);
          
          if ( this->FeatureEdgeSmoothing &&
               vtkMath::Dot(normal,neiNormal) <= CosFeatureAngle ) 
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }
        else // a visited edge; skip rest of analysis
        {
          continue;
        }
        
        if ( edge && Verts[p1].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p1].edges->Reset();
          Verts[p1].edges->InsertNextId(p2);
          Verts[p1].type = edge;
        }
        else if ( (edge && Verts[p1].type == VTK_BOUNDARY_EDGE_VERTEX) ||
                  (edge && Verts[p1].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p1].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p1].edges->InsertNextId(p2);
          if ( Verts[p1].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p1].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
        
        if ( edge && Verts[p2].type == VTK_SIMPLE_VERTEX )
        {
          Verts[p2].edges->Reset();
          Verts[p2].edges->InsertNextId(p1);
          Verts[p2].type = edge;
        }
        else if ( (edge && Verts[p2].type == VTK_BOUNDARY_EDGE_VERTEX ) ||
                  (edge && Verts[p2].type == VTK_FEATURE_EDGE_VERTEX) ||
                  (!edge && Verts[p2].type == VTK_SIMPLE_VERTEX ) )
        {
          Verts[p2].edges->InsertNextId(p1);
          if ( Verts[p2].type && edge == VTK_BOUNDARY_EDGE_VERTEX )
          {
            Verts[p2].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
      }
    }
    
    inMesh->Delete();
    if (toTris) {toTris->Delete();}
    
    neighbors->Delete();
  }//if strips or polys
  
    //post-process edge vertices to make sure we can smooth them
  for (i=0; i<numPts; i++)
  {
    if ( Verts[i].type == VTK_SIMPLE_VERTEX )
    {
      numSimple++;
    }
    
    else if ( Verts[i].type == VTK_FIXED_VERTEX )
    {
      numFixed++;
    }
    
    else if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX ||
              Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
    { //see how many edges; if two, what the angle is
      
      if ( !this->BoundarySmoothing && 
           Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX )
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numBEdges++;
      }
      
      else if ( (npts = Verts[i].edges->GetNumberOfIds()) != 2 )
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numFixed++;
      }
      
      else //check angle between edges
      {
        inPts->GetPoint(Verts[i].edges->GetId(0),x1);
        inPts->GetPoint(i,x2);
        inPts->GetPoint(Verts[i].edges->GetId(1),x3);
        
        for (k=0; k<3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize(l1) >= 0.0 &&
             vtkMath::Normalize(l2) >= 0.0 &&
             vtkMath::Dot(l1,l2) < CosEdgeAngle)
        {
          numFixed++;
          Verts[i].type = VTK_FIXED_VERTEX;
        }
        else
        {
          if ( Verts[i].type == VTK_FEATURE_EDGE_VERTEX )
          {
            numFEdges++;
          }
          else
          {
            numBEdges++;
          }
        }
      }//if along edge
    }//if edge vertex
  }//for all points
  
  if(DebugLevel>5) {
    cout<<"Found\n\t" << numSimple << " simple vertices\n\t"
      << numFEdges << " feature edge vertices\n\t"
      << numBEdges << " boundary edge vertices\n\t"
      << numFixed << " fixed vertices\n\t"<<endl;
    cout<<"1:numPts="<<numPts<<endl;
  }
  
  for (i=0; i<numPts; i++) 
  {
    if(DebugLevel>5) cout<<"Verts["<<i<<"].type="<<VertexType2Str(Verts[i].type)<<endl;
    if(Verts[i].edges != NULL && (npts = Verts[i].edges->GetNumberOfIds()) > 0)
    {
      for (j=0; j<npts; j++)
      {
        if(DebugLevel>5) cout<<"Verts["<<i<<"].edges->GetId("<<j<<")="<<Verts[i].edges->GetId(j)<<endl;
      }//for all connected points
    }
  }
  
  //Copy node type info from Verts
  EG_VTKDCN(vtkIntArray, node_type, m_grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    cout<<"Verts["<<node<<"].type="<<Verts[node].type<<endl;
    node_type->SetValue(node,Verts[node].type);
  }
  
  //free up connectivity storage
  for (int i=0; i<numPts; i++)
  {
    if ( Verts[i].edges != NULL )
    {
      Verts[i].edges->Delete();
      Verts[i].edges = NULL;
    }
  }
  delete [] Verts;
  
  return(0);
}
