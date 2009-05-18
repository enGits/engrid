#include "surfacemesher.h"

#include "insertpoints.h"
#include "removepoints.h"
#include "updatedesiredmeshdensity.h"

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
    
    //Method 3
    bool DEBUG=true;
    QString DEBUGDIR="/data1/home/mtaverne/Geometries/DEBUG/";
    
    if(insert_FP) {
      MeshDensityFunction();
      _insert_FP_all();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-swap-1");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-laplace");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-swap-2");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_FP");
    
    if(insert_EP) {
      MeshDensityFunction();
      _insert_EP_all();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-laplace");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_EP");
    
/*    if(insert_FP) {
      MeshDensityFunction();
      InsertPoints insert_field_points;
      insert_field_points.SetBCS(m_bcs);
      insert_field_points.Set_insert_FP(true);
      insert_field_points.Set_insert_EP(false);
      insert_field_points();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-swap-1");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-laplace");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-swap-2");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_FP");
    
    if(insert_EP) {
      MeshDensityFunction();
      InsertPoints insert_edge_points;
      insert_edge_points.SetBCS(m_bcs);
      insert_edge_points.Set_insert_FP(false);
      insert_edge_points.Set_insert_EP(true);
      insert_edge_points();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_EP-post-laplace");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_EP");*/
    
    if(remove_FP) {
      MeshDensityFunction();
      RemovePoints remove_field_points;
      remove_field_points.SetBCS(m_bcs);
      remove_field_points.Set_remove_FP(true);
      remove_field_points.Set_remove_EP(false);
      remove_field_points();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-remove");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-swap-1");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-laplace");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-swap-2");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-remove_FP");
    
    if(remove_EP) {
      MeshDensityFunction();
      RemovePoints remove_edge_points;
      remove_edge_points.SetBCS(m_bcs);
      remove_edge_points.Set_remove_FP(false);
      remove_edge_points.Set_remove_EP(true);
      remove_edge_points();
      if(DEBUG) DualSave(DEBUGDIR+"remove_EP-post-remove");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_EP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_EP-post-laplace");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-remove_EP");
    
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
  
  MeshDensityFunction();
  
  UpdateMeshDensity();
  if(i_iter<NumberOfIterations) cout<<"WARNING: Exited before finishing all iterations."<<endl;
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  
}
//end of operate()

void SurfaceMesher::MeshDensityFunction()
{
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setMaxiterDensity(MaxiterDensity);
  update_desired_mesh_density.SetVertexMeshDensityVector(VMDvector);
  update_desired_mesh_density();
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

//========================================================

bool SurfaceMesher::_insert_fieldpoint(vtkIdType D)
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

bool SurfaceMesher::_insert_edgepoint(vtkIdType j,vtkIdType K)// node1 K, node2 j
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

int SurfaceMesher::_insert_FP_counter()
{
  cout<<"===insert_FP_counter() START==="<<endl;
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
    if( !m_marked_cells[id_cell] && _insert_fieldpoint(id_cell) )
    {
      if(DebugLevel>0) cout<<"inserting a field point "<<id_cell<<endl;
      N_inserted_FP++;
      m_marked_cells[id_cell]=true;
      N_newcells+=2;
      N_newpoints+=1;
    }
  }
  cout<<"===insert_FP_counter() END==="<<endl;
  return(0);
}

int SurfaceMesher::_insert_EP_counter()
{
  cout<<"===insert_EP_counter() START==="<<endl;
  
  //Phase C: Prepare m_edge_map
  cout<<"===Phase C==="<<endl;
  m_edge_map.clear();
  vtkIdType edgeId=1;
  foreach(vtkIdType node1,m_SelectedNodes)
  {
//       cout<<"node1="<<node1<<endl;
    foreach(vtkIdType node2,n2n[node1])
    {
      if(m_edge_map[OrderedPair(node1,node2)]==0) { //this edge hasn't been numbered yet
        m_edge_map[OrderedPair(node1,node2)]=edgeId;edgeId++;
      }
    }
  }
  cout<<"m_edge_map.size()="<<m_edge_map.size()<<endl;
  
  m_StencilVector.clear();
  QMapIterator< pair<vtkIdType,vtkIdType>, vtkIdType> m_edge_map_iter(m_edge_map);
      //rewind the iterator
  m_edge_map_iter.toFront ();
      //start loop
  while (m_edge_map_iter.hasNext()) {
    m_edge_map_iter.next();
    vtkIdType node1=m_edge_map_iter.key().first;
    vtkIdType node2=m_edge_map_iter.key().second;
    if(DebugLevel>10) cout << "--->(" << node1 << "," << node2 << ")" << ": " << m_edge_map_iter.value() << endl;
    QSet <int> stencil_cells_set;
    QVector <int> stencil_cells_vector;
    stencil_cells_set=n2c[node1];
    stencil_cells_set.intersect(n2c[node2]);
    if(DebugLevel>10) cout<<"stencil_cells_set="<<stencil_cells_set<<endl;
    
    stencil_cells_vector.resize(stencil_cells_set.size());
    qCopy(stencil_cells_set.begin(),stencil_cells_set.end(),stencil_cells_vector.begin());
    if(DebugLevel>10) cout<<"stencil_cells_vector="<<stencil_cells_vector<<endl;
    
    vtkIdType id_cell=stencil_cells_vector[0];
    int SideToSplit = getSide(id_cell,grid,node1,node2);
    if(DebugLevel>10) cout<<"SideToSplit="<<SideToSplit<<endl;
    if(DebugLevel>10) cout<<"c2c[id_cell][SideToSplit]="<<c2c[id_cell][SideToSplit]<<endl;
    if(DebugLevel>10) for(int i=0;i<3;i++) cout<<"c2c[id_cell]["<<i<<"]="<<c2c[id_cell][i]<<endl;
    stencil_t S=getStencil(id_cell,SideToSplit);
    
    bool stencil_marked=false;
    foreach(vtkIdType C,stencil_cells_vector)
    {
      if(m_marked_cells[C]) stencil_marked=true;
    }
    if(DebugLevel>10) cout<<"stencil_marked="<<stencil_marked<<endl;
    if(DebugLevel>10) cout<<"insert_edgepoint(node1,node2)="<<_insert_edgepoint(node1,node2)<<endl;
    
    if( !stencil_marked && _insert_edgepoint(node1,node2) )
    {
      if(DebugLevel>1) cout<<"inserting an edge point "<< "(" << node1 << "," << node2 << ")" << ": " << m_edge_map_iter.value() << endl;
      N_inserted_EP++;
      foreach(vtkIdType C,stencil_cells_vector) m_marked_cells[C]=true;
      m_StencilVector.push_back(S);
      
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

int SurfaceMesher::_insert_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_FP_actor START==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
/*    if(m_marked_cells[id_cell]) cout<<"--->m_marked_cells["<<id_cell<<"]=TRUE"<<endl;
    else cout<<"--->m_marked_cells["<<id_cell<<"]=FALSE"<<endl;*/
    
    if( !m_marked_cells[id_cell] && _insert_fieldpoint(id_cell) )
    {
      if(DebugLevel>0) cout<<"inserting a field point "<<id_cell<<endl;
      vtkIdType newBC=cell_code_tmp->GetValue(id_cell);
      if(DebugLevel>42) cout<<"id_cell="<<id_cell<<" newBC="<<newBC<<endl;
      
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t C(0,0,0);
      
      int N_neighbours=N_pts;
      if(DebugLevel>42) cout<<"N_neighbours="<<N_neighbours<<endl;
      vec3_t corner[4];
      vtkIdType pts_triangle[4][3];
      for(int i=0;i<N_neighbours;i++)
      {
        grid->GetPoints()->GetPoint(pts[i], corner[i].data());
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

int SurfaceMesher::_insert_EP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_EP_actor START==="<<endl;
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(stencil_t S,m_StencilVector)
  {
    if(DebugLevel>10) cout<<"S="<<S<<endl;
    vec3_t A,B;
    grid_tmp->GetPoint(S.p[1],A.data());
    grid_tmp->GetPoint(S.p[3],B.data());
    vec3_t M=0.5*(A+B);
    addPoint(grid_tmp,m_newNodeId,M.data());
    if(DebugLevel>0) cout<<"NEW EDGE POINT: "<<m_newNodeId<<endl;
    
    vtkIdType pts_triangle[4][3];
    
    if(S.valid){//there is a neighbour cell
      if(DebugLevel>10) cout<<"m_marked_cells["<<S.id_cell1<<"]=true;"<<endl;
      if(DebugLevel>10) cout<<"m_marked_cells["<<S.id_cell2<<"]=true;"<<endl;
      m_marked_cells[S.id_cell1]=true;
      m_marked_cells[S.id_cell2]=true;
      
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
      if(DebugLevel>10) cout<<"m_marked_cells["<<S.id_cell1<<"]=true;"<<endl;
      m_marked_cells[S.id_cell1]=true;
      
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

int SurfaceMesher::_insert_FP_all()
{
  cout<<"===insert_FP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_inserted_FP=0;
  
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  _insert_FP_counter();
  
    //unmark cells (TODO: optimize)
  m_marked_cells.clear();
    //init grid_tmp
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  _insert_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,grid);
  cout<<"===insert_FP_all END==="<<endl;
  return(0);
}

int SurfaceMesher::_insert_EP_all()
{
  cout<<"===insert_EP_all START==="<<endl;
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity, grid, "node_meshdensity");
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  N_inserted_EP=0;
  
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  N_newpoints=0;
  N_newcells=0;
  
  m_marked_cells.clear();
  m_marked_nodes.clear();
  
  _insert_EP_counter();
  
    //unmark cells (TODO: optimize)
  m_marked_cells.clear();
    //init grid_tmp
  N_points=grid->GetNumberOfPoints();
  N_cells=grid->GetNumberOfCells();
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  m_total_N_newpoints+=N_newpoints; m_total_N_newcells+=N_newcells;
  
  makeCopyNoAlloc(grid, grid_tmp);
    //initialize new node counter
  m_newNodeId=N_points;
  
  _insert_EP_actor(grid_tmp);
  
  makeCopy(grid_tmp,grid);
  
  cout<<"===insert_EP_all END==="<<endl;
  return(0);
}
