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

bool RemovePoints::removePointCriteria(vtkIdType id_node)
{
  double QL1max=0.8;
  double QL2max=0.5;
  bool result1 = Q_L1(id_node)<QL1max && Q_L2(id_node)<QL2max;
  
  QVector <vtkIdType> PSP = getPotentialSnapPoints(id_node);
  double Lmean = CurrentVertexAvgDist(id_node);
  bool result2 = Lmean<desiredEdgeLength(PSP[0]) && Lmean<desiredEdgeLength(PSP[1]);
  
  return ( result1 );
}

void RemovePoints::operate()
{
  int N1 = grid->GetNumberOfPoints();
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  getNodesFromCells(m_SelectedCells, m_SelectedNodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  l2l_t n2c   = getPartN2C();
  l2l_t n2n   = getPartN2N();
  UpdateNodeType();
  
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
  
  int l_N_removed_FP = 0;
  
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  foreach(vtkIdType node,m_SelectedNodes)
  {
    if(node_type->GetValue(node)!=VTK_FIXED_VERTEX) {
      bool marked=false;
      foreach(vtkIdType id_cell,n2c[node]) {
        if(m_marked_cells[id_cell]) marked=true;
      }
      
      QSet <vtkIdType> DeadCells;
      QSet <vtkIdType> MutatedCells;
      QSet <vtkIdType> MutilatedCells;
      if( !marked && removePointCriteria(node) && FindSnapPoint(node,DeadCells,MutatedCells,MutilatedCells, N_newpoints, N_newcells)!=-1)
      {
        l_N_removed_FP++;
        m_hitlist[node]=1;
        foreach(vtkIdType id_cell,n2c[node]) m_marked_cells[id_cell]=true;
      }
    }
  }
  
  QSet <vtkIdType> DeadNodes;
  for(vtkIdType i=0;i<m_hitlist.size();i++) {
    if(m_hitlist[i]==2) DeadNodes.insert(i);
  }
  int N_newpoints=0;
  int N_newcells=0;
  DeleteSetOfPoints(DeadNodes, N_newpoints, N_newcells);
  
  int kills=-N_newpoints;
  int contracts=DeadNodes.size();
  if (kills != contracts) {
    EG_BUG;
  }
  
  int N2 = grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
}
