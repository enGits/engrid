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

SurfaceMesher::SurfaceMesher() : SurfaceOperation()
{
  EG_TYPENAME;
}

void SurfaceMesher::computeMeshDensity()
{
  ///@@@  TODO: Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(grid);
  update_desired_mesh_density.setConvergence_meshdensity(1e-4);
  update_desired_mesh_density.setMaxiterDensity(1000);
  update_desired_mesh_density.setVertexMeshDensityVector(VMDvector);
  update_desired_mesh_density();

  //UpdateCurrentMeshDensity();
  //UpdateNodeType();
}

void SurfaceMesher::updateNodeInfo(bool update_type)
{
  l2g_t nodes = getPartNodes();
  setAllCells();
  foreach(vtkIdType node, nodes) {
    if(update_type) {
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
}

void SurfaceMesher::swap()
{
  SwapTriangles swap;
  swap.setGrid(grid);
  swap.setRespectBC(true);
  swap.setFeatureSwap(true);
  QSet<int> rest_bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(rest_bcs);
  rest_bcs -= m_BCs;
  swap.setBoundaryCodes(rest_bcs);
  swap();
}

void SurfaceMesher::smooth(int N_iter)
{
  LaplaceSmoother lap;
  lap.setGrid(grid);
  QVector<vtkIdType> cls;
  getSurfaceCells(m_BCs, cls, grid);
  lap.setCells(cls);
  lap.setNumberOfIterations(N_iter);
  lap();
}

int SurfaceMesher::insertNodes()
{
  InsertPoints insert_points;
  insert_points.setGrid(grid);
  insert_points.setBCS(m_BCs);
  insert_points.set_insert_FP(false);
  insert_points.set_insert_EP(true);
  insert_points.setVertexMeshDensityVector(VMDvector);
  insert_points();
  return insert_points.getNumInserted();
}

int SurfaceMesher::deleteNodes()
{
  RemovePoints remove_points;
  remove_points.setGrid(grid);
  remove_points.setBCS(m_BCs);
  remove_points.set_remove_FP(true);
  remove_points.set_remove_EP(true);
  remove_points();
  return remove_points.getNumRemoved();
}

void SurfaceMesher::operate()
{
  updateNodeInfo(true);
  int num_inserted = 0;
  int num_deleted = 0;
  do {
    computeMeshDensity();
    num_inserted = insertNodes();
    //cout << num_inserted << " nodes inserted" << endl;
    swap();
    smooth(4);
    /*
    do {
      num_deleted = deleteNodes();
      cout << num_deleted << " nodes deleted" << endl;
    } while (num_deleted > 0);
    */
    for (int i = 0; i < 10; ++i) {
      swap();
      smooth(1);
    }
  } while (num_inserted - num_deleted > grid->GetNumberOfPoints()/100);
  //createIndices(grid);
  updateNodeInfo(true);
}

