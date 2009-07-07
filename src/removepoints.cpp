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
  setQuickSave( true );
}

bool RemovePoints::removePointCriteria( vtkIdType id_node )
{
  double QL1max = 0.8;
  double QL2max = 0.5;
  bool result1 = Q_L1( id_node ) < QL1max && Q_L2( id_node ) < QL2max;

  QVector <vtkIdType> PSP = getPotentialSnapPoints( id_node );
  double Lmean = CurrentVertexAvgDist( id_node );
  bool result2 = Lmean < desiredEdgeLength( PSP[0] ) && Lmean < desiredEdgeLength( PSP[1] );

  return ( result1 );
}

void RemovePoints::operate()
{
  int N1 = grid->GetNumberOfPoints();

  QVector<vtkIdType> selected_cells;
  getSurfaceCells( m_bcs, selected_cells, grid );
  QVector <vtkIdType> selected_nodes;
  getNodesFromCells( selected_cells, selected_nodes, grid );

  setAllSurfaceCells();
  l2l_t n2c   = getPartN2C();

  UpdatePotentialSnapPoints(true);

  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );

  //for FindSnapPoint + DeleteSetOfPoints
  int num_newpoints = 0;
  int num_newcells = 0;
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;

  QMap <vtkIdType, bool> marked_cells;
  QMap <vtkIdType, bool> marked_nodes;
  QSet <vtkIdType> DeadNodes;

  //count
  foreach( vtkIdType node, selected_nodes ) {
    if ( node_type->GetValue( node ) != VTK_FIXED_VERTEX ) {
      bool marked = false;
      foreach( vtkIdType id_cell, n2c[node] ) {
        if ( marked_cells[id_cell] ) marked = true;
      }

      if ( !marked && removePointCriteria( node ) && FindSnapPoint( node, DeadCells, MutatedCells, MutilatedCells, num_newpoints, num_newcells ) != -1 ) {
        DeadNodes.insert( node );
        foreach( vtkIdType id_cell, n2c[node] ) marked_cells[id_cell] = true;
      }
    }
  }

  //delete
  DeleteSetOfPoints( DeadNodes, num_newpoints, num_newcells );

  int N2 = grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
}
