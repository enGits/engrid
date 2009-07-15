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
  getSet("surface meshing", "maximal number of iterations", 20, m_NumMaxIter);
  getSet("surface meshing", "number of smoothing steps", 1, m_NumSmoothSteps);
}

void SurfaceMesher::computeMeshDensity()
{
  ///@@@  TODO: Optimize by using only one loop through nodes!
  UpdateDesiredMeshDensity update_desired_mesh_density;
  update_desired_mesh_density.setGrid(grid);
  update_desired_mesh_density.setVertexMeshDensityVector(VMDvector);
  update_desired_mesh_density();
}

void SurfaceMesher::updateNodeInfo(bool update_type)
{
  setAllCells();
  l2g_t nodes = getPartNodes();
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
      //node_meshdensity_desired->SetValue(node, VMDvector[idx].density);
    } else { //unspecified
      //double D=DesiredMeshDensity(node);
      //node_meshdensity_desired->SetValue(node, D);
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
  insert_points.setBoundaryCodes(m_BCs);
  insert_points();
  return insert_points.getNumInserted();
}

int SurfaceMesher::deleteNodes()
{
  RemovePoints remove_points;
  remove_points.setGrid(grid);
  remove_points.setBCS(m_BCs);
  remove_points();
  return remove_points.getNumRemoved();
}

void SurfaceMesher::operate()
{
  EG_VTKDCN(vtkDoubleArray, md, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    md->SetValue(id_node, 1e-6);
  }
  updateNodeInfo(true);
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  while (!done) {
    computeMeshDensity();
    num_inserted = insertNodes();
    swap();
    num_deleted = 0;
    int N = 0;
    int count = 0;

    do {
      N = deleteNodes();
      num_deleted += N;
      ++count;
    } while ((N > 0) && (count < 20));

    for (int i = 0; i < m_NumSmoothSteps; ++i) {
      swap();
      smooth(1);
    }
    ++iter;
    done = (iter >= m_NumMaxIter) || (num_inserted - num_deleted < grid->GetNumberOfPoints()/100);
    cout << "surface mesher iteration " << iter << ":" << endl;
    cout << "  inserted nodes : " << num_inserted << endl;
    cout << "  deleted nodes  : " << num_deleted << endl;
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  createIndices(grid);
  updateNodeInfo(true);
  computeMeshDensity();
}
