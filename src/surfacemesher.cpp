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
    cout<<"surface meshing iteration " << i_iter << "/" << NumberOfIterations << endl;
    
    m_total_N_newpoints=0;
    m_total_N_newcells=0;
    
    setAllSurfaceCells();
    
    //edit points
    //cout<<"===Phase D==="<<endl;
    N_inserted_FP=0;
    N_inserted_EP=0;
    N_removed_FP=0;
    N_removed_EP=0;
    
    //Method 3
    if (insert_FP) {
      UpdateNodeInfo(false);
      InsertPoints insert_field_points;
      insert_field_points.setGrid(grid);
      insert_field_points.setBCS(m_bcs);
      insert_field_points.set_insert_FP(true);
      insert_field_points.set_insert_EP(false);
      insert_field_points.setVertexMeshDensityVector(VMDvector);
      insert_field_points();
      if (DoSwap) {
        SwapFunction();
      }
      if (DoLaplaceSmoothing) {
        SmoothFunction();
      }
      MeshDensityFunction();
    }
    
    if (insert_EP) {
//    MeshDensityFunction();
      UpdateNodeInfo(false);
      InsertPoints insert_edge_points;
      insert_edge_points.setGrid(grid);
      insert_edge_points.setBCS(m_bcs);
      insert_edge_points.set_insert_FP(false);
      insert_edge_points.set_insert_EP(true);
      insert_edge_points.setVertexMeshDensityVector(VMDvector);
      insert_edge_points();
      if (DoSwap) {
        SwapFunction();
      }
      MeshDensityFunction();
    }
    
    if (remove_FP) {
      int N;
      do {
        UpdateNodeInfo(false);
        RemovePoints remove_field_points;
        remove_field_points.setGrid(grid);
        remove_field_points.setBCS(m_bcs);
        remove_field_points.set_remove_FP(true);
        remove_field_points.set_remove_EP(false);
        remove_field_points();
        SwapFunction();
      } while (N > 0);
      MeshDensityFunction();
    }
    
    if (remove_EP) {
      int N;
      do {
        UpdateNodeInfo(false);
        RemovePoints remove_edge_points;
        remove_edge_points.setGrid(grid);
        remove_edge_points.setBCS(m_bcs);
        remove_edge_points.set_remove_FP(false);
        remove_edge_points.set_remove_EP(true);
        remove_edge_points();
        N = remove_edge_points.getNumRemoved();
        SwapFunction();
      } while (N > 0);
      MeshDensityFunction();
    }
    
    if (DoLaplaceSmoothing) {
      SmoothFunction();
    }
    MeshDensityFunction();

    cout << "  cells: " << grid->GetNumberOfCells() << endl;
    cout << "  nodes: " << grid->GetNumberOfPoints() << endl;

  }
  
  if (i_iter < NumberOfIterations) {
    cout << " WARNING: Exited before finishing all iterations." << endl;
  }
  
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

  //UpdateCurrentMeshDensity();
  //UpdateNodeType();
}

void SurfaceMesher::UpdateNodeInfo(bool UpdateType)
{
  l2g_t nodes = getPartNodes();

  //cout<<"=== UpdateNodeInfo START ==="<<endl;
  setAllCells();
  foreach(vtkIdType node, nodes) {
    if(UpdateType) {
      static int nStatic_UpdateType;    // Value of nStatic_UpdateType is retained between each function call
      nStatic_UpdateType++;
      //cout << "nStatic_UpdateType is " << nStatic_UpdateType << endl;
      
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
  //cout<<"=== UpdateNodeInfo STOP ==="<<endl;
}

int SurfaceMesher::SwapFunction()
{
  //cout<<"=== SwapFunction START ==="<<endl;
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
  
  //cout<<"=== SwapFunction END ==="<<endl;
  return(0);
}

int SurfaceMesher::SmoothFunction()
{
  LaplaceSmoother lap;
  lap.setGrid(this->grid);
  QVector<vtkIdType> cls;
  getSurfaceCells(m_bcs, cls, this->grid);
  lap.setCells(cls);
  lap.setNumberOfIterations(N_SmoothIterations);
  lap();
  return(0);
}
