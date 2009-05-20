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
//       MeshDensityFunction();
      InsertPoints insert_field_points;
      insert_field_points.SetBCS(m_bcs);
      insert_field_points.Set_insert_FP(true);
      insert_field_points.Set_insert_EP(false);
      insert_field_points();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-insert");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"insert_FP-post-laplace");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_FP");
    
    if(insert_EP) {
//       MeshDensityFunction();
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
    if(DEBUG) DualSave(DEBUGDIR+"post-insert_EP");
    
    if(remove_FP) {
//       MeshDensityFunction();
      RemovePoints remove_field_points;
      remove_field_points.SetBCS(m_bcs);
      remove_field_points.Set_remove_FP(true);
      remove_field_points.Set_remove_EP(false);
      remove_field_points();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-remove");
      if(DoSwap) SwapFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-swap");
      if(DoLaplaceSmoothing) SmoothFunction();
      if(DEBUG) DualSave(DEBUGDIR+"remove_FP-post-laplace");
    }
    if(DEBUG) DualSave(DEBUGDIR+"post-remove_FP");
    
    if(remove_EP) {
//       MeshDensityFunction();
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
