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
#include "surfacemesher.h"

#include "insertpoints.h"
#include "removepoints.h"
#include "updatedesiredmeshdensity.h"

#include <vtkSmoothPolyDataFilter.h>

SurfaceMesher::SurfaceMesher()
: SurfaceOperation()
{
  EG_TYPENAME;
}

void SurfaceMesher::operate()
{
  QTime start = QTime::currentTime();
  
//   UpdateNodeInfo(true);
  MeshDensityFunction();
  
  int i_iter=0;
  ///@@@  TODO:Optimize this loop
  for(i_iter=0;i_iter<NumberOfIterations;i_iter++) {
    cout<<"===ITERATION NB "<<i_iter<<"/"<<NumberOfIterations<<"==="<<endl;
    
    m_total_N_newpoints=0;
    m_total_N_newcells=0;
    
    setAllSurfaceCells();
    
    //edit points
    cout<<"===Phase D==="<<endl;
    N_inserted_FP=0;
    N_inserted_EP=0;
    N_removed_FP=0;
    N_removed_EP=0;
    
    //Method 3
    if(insert_FP) {
//       MeshDensityFunction();
      UpdateNodeInfo(false);
      InsertPoints insert_field_points;
      insert_field_points.setGrid(grid);
      insert_field_points.setSource(m_ProjectionSurface);
      insert_field_points.setBCS(m_bcs);
      insert_field_points.set_insert_FP(true);
      insert_field_points.set_insert_EP(false);
      insert_field_points.setVertexMeshDensityVector(VMDvector);
      insert_field_points();
      insert_field_points.delete_CellLocator_and_ProjectionSurface();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(insert_EP) {
//       MeshDensityFunction();
      UpdateNodeInfo(false);
      InsertPoints insert_edge_points;
      insert_edge_points.setGrid(grid);
      insert_edge_points.setSource(m_ProjectionSurface);
      insert_edge_points.setBCS(m_bcs);
      insert_edge_points.set_insert_FP(false);
      insert_edge_points.set_insert_EP(true);
      insert_edge_points.setVertexMeshDensityVector(VMDvector);
      insert_edge_points();
      insert_edge_points.delete_CellLocator_and_ProjectionSurface();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(remove_FP) {
//       MeshDensityFunction();
      UpdateNodeInfo(false);
      RemovePoints remove_field_points;
      remove_field_points.setGrid(grid);
      remove_field_points.setBCS(m_bcs);
      remove_field_points.set_remove_FP(true);
      remove_field_points.set_remove_EP(false);
      remove_field_points();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
    if(remove_EP) {
//       MeshDensityFunction();
      UpdateNodeInfo(false);
      RemovePoints remove_edge_points;
      remove_edge_points.setGrid(grid);
      remove_edge_points.setBCS(m_bcs);
      remove_edge_points.set_remove_FP(false);
      remove_edge_points.set_remove_EP(true);
      remove_edge_points();
      if(DoSwap) SwapFunction();
      if(DoLaplaceSmoothing) SmoothFunction();
    }
    
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
//     if(N_inserted_FP==0 && N_inserted_EP==0 && N_removed_FP==0 && N_removed_EP==0) break;

    cout << "i_iter/NumberOfIterations=" << i_iter << "/" << NumberOfIterations << endl;
  }
  
  cout << "i_iter/NumberOfIterations=" << i_iter << "/" << NumberOfIterations << endl;
  
/*  MeshDensityFunction();
  UpdateCurrentMeshDensity();*/
  
  if(i_iter<NumberOfIterations) cout<<"WARNING: Exited before finishing all iterations."<<endl;
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  
}
//end of operate()

void SurfaceMesher::MeshDensityFunction()
{
  ///@@@  TODO: Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(grid);
  update_desired_mesh_density.setConvergence_meshdensity(Convergence_meshdensity);
  update_desired_mesh_density.setMaxiterDensity(MaxiterDensity);
  update_desired_mesh_density.setVertexMeshDensityVector(VMDvector);
  update_desired_mesh_density();
/*  UpdateCurrentMeshDensity();
  UpdateNodeType();*/
}

void SurfaceMesher::UpdateNodeInfo(bool UpdateType)
{
  l2g_t nodes = getPartNodes();

  cout<<"=== UpdateNodeInfo START ==="<<endl;
  setAllCells();
  foreach(vtkIdType node, nodes) {
    if(UpdateType) {
      static int nStatic_UpdateType;    // Value of nStatic_UpdateType is retained between each function call
      nStatic_UpdateType++;
      cout << "nStatic_UpdateType is " << nStatic_UpdateType << endl;
      
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      node_type->SetValue(node, getNodeType(node));
    }
    EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");//what we have
    node_meshdensity_current->SetValue(node, CurrentMeshDensity(node));
    
    EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");//density index from table
    VertexMeshDensity nodeVMD = getVMD(node);
    int idx=VMDvector.indexOf(nodeVMD);
    node_specified_density->SetValue(node, idx);
    
    EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");//what we want
    if(idx!=-1) { //specified
      node_meshdensity_desired->SetValue(node, VMDvector[idx].density);
    } else { //unspecified
      double D=DesiredMeshDensity(node);
      node_meshdensity_desired->SetValue(node, D);
    }
  }
  cout<<"=== UpdateNodeInfo STOP ==="<<endl;
}

int SurfaceMesher::SwapFunction()
{
  cout<<"=== SwapFunction START ==="<<endl;
  QSet<int> rest_bcs;
  mainWindow()->getAllBoundaryCodes(rest_bcs);
  rest_bcs -= m_bcs;
  GuiMainWindow::pointer()->quickSave(GuiMainWindow::pointer()->getFilePath() + "beforeswap");
  setAllSurfaceCells();
  SwapTriangles swap;
  swap.setGrid(this->grid);
  swap.setQuickSave(true);
  swap.setRespectBC(true);
  swap.setFeatureSwap(true);
  swap.setBoundaryCodes(rest_bcs);
  swap();
  
  cout<<"=== SwapFunction END ==="<<endl;
  return(0);
}

int SurfaceMesher::SmoothFunction()
{
  cout<<"=== SmoothFunction START ==="<<endl;
  UpdateNodeInfo(false);
  // translate points to smooth grid
  //4 possibilities
  //vtk smooth 1
  //vtk smooth 2
  //laplacian smoothing with projection
  //Roland smoothing with projection
  
  if (true) {
    //laplacian smoothing with projection
    LaplaceSmoother lap;
    lap.setGrid(this->grid);
    QVector<vtkIdType> cls;
    getSurfaceCells(m_bcs, cls, this->grid);
    lap.setCells(cls);
    lap.setNumberOfIterations(N_SmoothIterations);
    lap();
  } else {
    //preparations
    QVector<vtkIdType> cells;
    getSurfaceCells(m_bcs, cells, this->grid);
    EG_VTKSP(vtkPolyData, pdata);
    addToPolyData(cells, pdata, this->grid);
    EG_VTKSP(vtkSmoothPolyDataFilter, smooth);
    
    //configure vtkSmoothPolyDataFilter
    smooth->SetInput(pdata);
    
    // smooth->SetConvergence (ui.doubleSpinBox_Convergence->value());
    smooth->SetNumberOfIterations (N_SmoothIterations);
    // smooth->SetRelaxationFactor (ui.lineEdit_RelaxationFactor->text().toDouble());
    smooth->SetFeatureEdgeSmoothing (false);
    // smooth->SetFeatureAngle (ui.doubleSpinBox_FeatureAngle->value());
    // smooth->SetEdgeAngle (ui.doubleSpinBox_EdgeAngle->value());
    smooth->SetBoundarySmoothing (true);
    // smooth->SetGenerateErrorScalars (ui.checkBox_GenerateErrorScalars->checkState());
    // smooth->SetGenerateErrorVectors (ui.checkBox_GenerateErrorVectors->checkState());
    
    QVector<vtkIdType> cells_Source;
    getAllSurfaceCells(cells_Source, m_ProjectionSurface);
    EG_VTKSP(vtkPolyData, pdata_Source);
    addToPolyData(cells_Source, pdata_Source, m_ProjectionSurface);
    smooth->SetSource (pdata_Source);
    
    //smooth
    smooth->Update();
    
    //copy smoothed grid to main grid
    EG_VTKDCN(vtkLongArray_t, node_index, pdata, "node_index");
    for (vtkIdType i = 0; i < smooth->GetOutput()->GetNumberOfPoints(); ++i) {
      vec3_t x;
      smooth->GetOutput()->GetPoints()->GetPoint(i, x.data());
      vtkIdType nodeId = node_index->GetValue(i);
      this->grid->GetPoints()->SetPoint(nodeId, x.data());
    }
  }
  
  cout<<"=== SmoothFunction END ==="<<endl;
  return(0);
}
