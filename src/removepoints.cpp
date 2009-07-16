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

#include "geometrytools.h"
using namespace GeometryTools;

#include <cmath>
using namespace std;

#include <iostream>
using namespace std;

RemovePoints::RemovePoints() : SurfaceOperation()
{
  setQuickSave( true );
  getSet("surface meshing", "point removal threshold", 1, m_Threshold);
}

///@@@ TODO: Replace sets with vectors if possible
void RemovePoints::operate()
{
  int N1 = grid->GetNumberOfPoints();

  QVector<vtkIdType> selected_cells;
  getSurfaceCells(m_BoundaryCodes, selected_cells, grid);
  QVector<vtkIdType> selected_nodes;
  getNodesFromCells(selected_cells, selected_nodes, grid);

  setAllSurfaceCells();
  l2l_t  n2n   = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  nodes = getPartNodes();
  
  UpdatePotentialSnapPoints(false);
  
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type" );
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code" );
  EG_VTKDCN(vtkDoubleArray, cl, grid, "node_meshdensity_desired" );

  // global values
  QSet <vtkIdType> all_deadcells;
  QSet <vtkIdType> all_mutatedcells;
  QSet <vtkIdType> all_mutilatedcells;
  int num_newpoints = 0;
  int num_newcells = 0;
  
  QVector <bool> marked_nodes(nodes.size(), false);
  
  QVector <vtkIdType> deadnode_vector;
  QVector <vtkIdType> snappoint_vector;
  
  //count
  for (int i_selected_nodes = 0; i_selected_nodes < selected_nodes.size(); ++i_selected_nodes) {
    vtkIdType id_node = selected_nodes[i_selected_nodes];
    int i_nodes = _nodes[id_node];
    if (node_type->GetValue(id_node) != VTK_FIXED_VERTEX) {
      if (!marked_nodes[i_nodes]) {
        vec3_t xi;
        grid->GetPoint(id_node, xi.data());
        double cl_node = cl->GetValue(id_node);
        bool remove_node = true;
        for (int j = 0; j < n2n[i_nodes].size(); ++j) {
          vtkIdType id_neigh = nodes[n2n[i_nodes][j]];
          double cl_neigh = cl->GetValue(id_neigh);
          vec3_t xj;
          grid->GetPoint(id_neigh, xj.data());
          double L = (xi-xj).abs();
          if (L > 0.5*(cl_node+cl_neigh)/m_Threshold) {
            remove_node = false;
            break;
          }
        }
        if (remove_node) {
          // local values
          QSet <vtkIdType> dead_cells;
          QSet <vtkIdType> mutated_cells;
          QSet <vtkIdType> mutilated_cells;
          int l_num_newpoints = 0;
          int l_num_newcells = 0;
          vtkIdType snap_point = FindSnapPoint(id_node, dead_cells, mutated_cells, mutilated_cells, l_num_newpoints, l_num_newcells);
          if(snap_point >= 0) {
            // add deadnode/snappoint pair
            deadnode_vector.push_back(id_node);
            snappoint_vector.push_back(snap_point);
            // update global values
            num_newpoints += l_num_newpoints;
            num_newcells  += l_num_newcells;
            all_deadcells.unite(dead_cells);
            all_mutatedcells.unite(mutated_cells);
            all_mutilatedcells.unite(mutilated_cells);
            // mark neighbour nodes
            foreach(int i_node_neighbour, n2n[_nodes[id_node]]) {
              marked_nodes[nodes[i_node_neighbour]] = true;
            }
          }
        }
      }
    }
  }

  //delete
  DeleteSetOfPoints(deadnode_vector, snappoint_vector, all_deadcells, all_mutatedcells, all_mutilatedcells, num_newpoints, num_newcells);

  int N2 = grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
}

// DEFINITIONS:
// Normal cell: nothing has changed
// Dead cell: the cell does not exist anymore
// Mutated cell: the cell's form has changed
// Mutilated cell: the cell has less points than before

///@@@  TODO: Clean up this function
vtkIdType RemovePoints::FindSnapPoint( vtkIdType DeadNode, QSet <vtkIdType> & DeadCells, QSet <vtkIdType> & MutatedCells, QSet <vtkIdType> & MutilatedCells, int& num_newpoints, int& num_newcells )
{
  // preparations
  l2l_t n2c = getPartN2C();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t cells = getPartCells();
  
  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  if ( node_type->GetValue( DeadNode ) == VTK_FIXED_VERTEX ) {
    cout << "ERROR: unable to remove fixed vertex." << endl;
    EG_BUG;
    return( -1 );
  }
  
  vtkIdType SnapPoint = -1;
  
  QVector <vtkIdType> PSP_vector = getPotentialSnapPoints( DeadNode );
  foreach( vtkIdType PSP, PSP_vector ) {
    bool IsValidSnapPoint = true;
    
    // TEST 0: DeadNode, PSP and any common point must belong to a cell.
    
    // TEST 1: Number of common points must not exceed 2.
    bool IsTetra = true;
    if ( NumberOfCommonPoints( DeadNode, PSP, IsTetra ) > 2 ) { //common point check
      if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
      IsValidSnapPoint = false;
    }
    // TEST 2: DeadNode, PSP and common points must not form a tetrahedron.
    if ( IsTetra ) { //tetra check
      if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
      IsValidSnapPoint = false;
    }
    
    //count number of points and cells to remove + analyse cell transformations
    num_newpoints = -1;
    num_newcells = 0;
    DeadCells.clear();
    MutatedCells.clear();
    MutilatedCells.clear();
    foreach( int i_cell, n2c[_nodes[DeadNode]] ) { //loop through potentially dead cells
      vtkIdType id_cell = cells[i_cell];
      //get points around cell
      vtkIdType num_pts, *pts;
      grid->GetCellPoints( id_cell, num_pts, pts );
      
      if ( num_pts != 3 ) {
        cout << "ERROR: Non-triangle detected!" << endl;
        EG_BUG;
      }
      
      bool ContainsSnapPoint = false;
      bool invincible = false;
      for ( int i = 0; i < num_pts; ++i ) {
        if ( pts[i] == PSP ) {
          ContainsSnapPoint = true;
        }
        if ( pts[i] != DeadNode && pts[i] != PSP &&  n2c[_nodes[pts[i]]].size() <= 1 ) {
          invincible = true;
        }
      }
      if ( ContainsSnapPoint ) { // potential dead cell
        if ( invincible ) {
          // TEST 3: Check that empty lines aren't left behind when a cell is killed
          if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
          IsValidSnapPoint = false;
        }
        else {
          DeadCells.insert( id_cell );
          num_newcells -= 1;
        }
      }
      else { // if the cell does not contain the SnapPoint (potential mutated cell)
        
        vtkIdType OldTriangle[3];
        vtkIdType NewTriangle[3];
        
        for ( int i = 0; i < num_pts; ++i ) {
          OldTriangle[i] = pts[i];
          NewTriangle[i] = (( pts[i] == DeadNode ) ? PSP : pts[i] );
        }
        vec3_t Old_N = triNormal( grid, OldTriangle[0], OldTriangle[1], OldTriangle[2] );
        vec3_t New_N = triNormal( grid, NewTriangle[0], NewTriangle[1], NewTriangle[2] );
        
        // TEST 4: area + inversion check
        if ( Old_N*New_N < 0 || New_N*New_N < Old_N*Old_N*1. / 100. ) {
          if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
          IsValidSnapPoint = false;
        }
        
        // TEST 5: flipped cell test from old laplace smoother
        vec3_t P;
        grid->GetPoint( PSP, P.data() );
        if ( FlippedCells( DeadNode, P ) ) {
          if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
          IsValidSnapPoint = false;
        }
        
        //mutated cell
        MutatedCells.insert( id_cell );
      }
    }
    
    // TEST 5: survivor check
    if ( grid->GetNumberOfCells() + num_newcells <= 0 ) {
      if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
      IsValidSnapPoint = false;
    }
    
    if ( IsValidSnapPoint ) {
      SnapPoint = PSP;
      break;
    }
  } //end of loop through potential SnapPoints
  
  return ( SnapPoint );
}
//End of FindSnapPoint

bool RemovePoints::DeleteSetOfPoints( QVector <vtkIdType> deadnode_vector, QVector <vtkIdType> snappoint_vector, QSet <vtkIdType> & all_deadcells, QSet <vtkIdType> & all_mutatedcells, QSet <vtkIdType> & all_mutilatedcells, int& num_newpoints, int& num_newcells)
{
  QVector <vtkIdType> inter_vector;
  QSet <vtkIdType> inter_set;
  qcontIntersection(deadnode_vector, snappoint_vector, inter_vector); if(inter_vector.size()>0) EG_BUG;
  qcontIntersection(all_deadcells, all_mutatedcells, inter_set); if(inter_set.size()>0) EG_BUG;
  qcontIntersection(all_deadcells, all_mutilatedcells, inter_set); if(inter_set.size()>0) EG_BUG;
  qcontIntersection(all_mutatedcells, all_mutilatedcells, inter_set); if(inter_set.size()>0) EG_BUG;
  
  int initial_num_points = grid->GetNumberOfPoints();
  
  //src grid info
  int num_points = grid->GetNumberOfPoints();
  int num_cells = grid->GetNumberOfCells();
  
  if ( num_newcells != 2*num_newpoints ) {
    EG_BUG;
  }
  
  //allocate
  EG_VTKSP( vtkUnstructuredGrid, dst );
  allocateGrid( dst, num_cells + num_newcells, num_points + num_newpoints );
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet( num_points );
  
  //copy undead points
  vtkIdType dst_id_node = 0;
  for ( vtkIdType src_id_node = 0; src_id_node < num_points; src_id_node++ ) {//loop through src points
    OffSet[src_id_node] = src_id_node - dst_id_node;
    if ( !deadnode_vector.contains( src_id_node ) ) { //if the node isn't dead, copy it
      vec3_t x;
      grid->GetPoints()->GetPoint( src_id_node, x.data() );
      dst->GetPoints()->SetPoint( dst_id_node, x.data() );
      copyNodeData( grid, src_id_node, dst, dst_id_node );
      dst_id_node++;
    }
    else {
      if ( DebugLevel > 0 ) {
        cout << "dead node encountered: src_id_node=" << src_id_node << " dst_id_node=" << dst_id_node << endl;
      }
    }
  }
  //Copy undead cells
  for ( vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell ) {//loop through src cells
    if ( !all_deadcells.contains( id_cell ) ) { //if the cell isn't dead
      vtkIdType src_num_pts, *src_pts;
      vtkIdType dst_num_pts, dst_pts[3];
      grid->GetCellPoints( id_cell, src_num_pts, src_pts );
      vtkIdType type_cell = grid->GetCellType( id_cell );
      
      dst_num_pts = 3;//src_num_pts;
      
      if ( all_mutatedcells.contains( id_cell ) ) { //mutated cell
        int num_deadnode = 0;
        for ( int i = 0; i < src_num_pts; i++ ) {
          int DeadIndex = deadnode_vector.indexOf( src_pts[i] );
          if ( DeadIndex != -1 ) {
            dst_pts[i] = snappoint_vector[DeadIndex] - OffSet[snappoint_vector[DeadIndex]]; // dead node
            num_deadnode++;
          }
          else {
            dst_pts[i] = src_pts[i] - OffSet[src_pts[i]]; // not a dead node
          }
        }
        if ( num_deadnode != 1 ) {
          qWarning() << "FATAL ERROR: Mutated cell has more than one dead node!";
          EG_BUG;
        }
      }
      else if ( all_mutilatedcells.contains( id_cell ) ) { //mutilated cell (ex: square becoming triangle) (WARNING: Not fully functional yet)
        cout << "FATAL ERROR: Quads not supported yet." << endl;
        EG_BUG;
        if ( type_cell == VTK_QUAD ) {
          type_cell = VTK_TRIANGLE;
          dst_num_pts -= 1;
        }
        else {
          cout << "FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!" << endl;
          EG_BUG;
        }
        // merge points
        for ( int i = 0; i < src_num_pts; i++ ) {
          ///@@@ TODO: finish this eventually for quad support
          //do nothing in case of deadnode_vector[i]
        }
      }
      else { //normal cell
        if ( DebugLevel > 10 ) {
          cout << "processing normal cell " << id_cell << endl;
        }
        for ( int i = 0; i < src_num_pts; i++ ) {
          int DeadIndex = deadnode_vector.indexOf( src_pts[i] );
          if ( DeadIndex != -1 ) {
            qWarning() << "FATAL ERROR: Normal cell contains a dead node!";
            EG_BUG;
          }
          dst_pts[i] = src_pts[i] - OffSet[src_pts[i]];
        }
      }
      // copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell( type_cell, dst_num_pts, dst_pts );
      copyCellData( grid, id_cell, dst, id_new_cell );
    }
  }
  
  makeCopy( dst, grid );
  
  if ( -num_newpoints != deadnode_vector.size() ) {
    EG_BUG;
  }
  
  int final_num_points = grid->GetNumberOfPoints();
  if ( initial_num_points - final_num_points != deadnode_vector.size() ) {
    EG_BUG;
  }
  
  return( true );
}
//End of DeleteSetOfPoints

int RemovePoints::NumberOfCommonPoints( vtkIdType id_node1, vtkIdType id_node2, bool& IsTetra )
{
  l2l_t  n2n   = getPartN2N();
  l2l_t  n2c   = getPartN2C();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes  = getPartNodes();
  l2g_t cells = getPartCells();
  
  QVector<int> node1_neighbours = n2n[_nodes[id_node1]];
  QVector<int> node2_neighbours = n2n[_nodes[id_node2]];
  QVector<int> intersection;
  qcontIntersection( node1_neighbours, node2_neighbours, intersection );
  int N = intersection.size();
  IsTetra = false;
  if ( N == 2 ) {
    vtkIdType intersection1 = nodes[intersection[0]];
    vtkIdType intersection2 = nodes[intersection[1]];
    
    // test if id_node1, id_node2 and intersection* form a cell
    QVector <vtkIdType> EdgeCells_1i;
    QVector <vtkIdType> EdgeCells_2i;
    QVector <vtkIdType> inter;
    int N;
    
    // intersection1
    N = getEdgeCells( id_node1, intersection1, EdgeCells_1i );
    if ( N != 2 ) EG_BUG;
    N = getEdgeCells( id_node2, intersection1, EdgeCells_2i );
    if ( N != 2 ) EG_BUG;
    qcontIntersection( EdgeCells_1i, EdgeCells_2i, inter );
    if ( inter.size() <= 0 ) EG_BUG;
    
    // intersection2
    N = getEdgeCells( id_node1, intersection2, EdgeCells_1i );
    if ( N != 2 ) EG_BUG;
    N = getEdgeCells( id_node2, intersection2, EdgeCells_2i );
    if ( N != 2 ) EG_BUG;
    qcontIntersection( EdgeCells_1i, EdgeCells_2i, inter );
    if ( inter.size() <= 0 ) EG_BUG;
    
    // check if DeadNode, PSP and common points form a tetrahedron.
    if ( n2n[_nodes[intersection1]].contains( _nodes[intersection2] ) ) { //if there's an edge between intersection1 and intersection2
      //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
      QVector<int> S1 = n2c[_nodes[intersection1]];
      QVector<int> S2 = n2c[_nodes[intersection2]];
      QVector<int> Si;
      qcontIntersection( S1, S2, Si );
      int counter = 0;
      foreach( int i_cell, Si ) {
        vtkIdType num_pts, *pts;
        grid->GetCellPoints( cells[i_cell], num_pts, pts );
        for ( int i = 0; i < num_pts; ++i ) {
          if ( pts[i] == id_node1 || pts[i] == id_node2 ) counter++;
        }
      }
      if ( counter >= 2 ) {
        IsTetra = true;
      }
    }
  }
  return( N );
}

bool RemovePoints::FlippedCells( vtkIdType id_node, vec3_t P )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();
  
  vec3_t x0_old, x0_new;
  grid->GetPoint( id_node, x0_old.data() );
  x0_new = P;
  
  foreach( int i_cell, n2c[_nodes[id_node]] ) {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType num_pts, *pts;
    grid->GetCellPoints( id_cell, num_pts, pts );
    int i;
    for ( i = 0; i < num_pts; i++ ) {
      if ( pts[i] == id_node ) {
        break;
      }
    }
    vec3_t x2, x3;
    grid->GetPoint( pts[( i+1 )%num_pts], x2.data() );
    grid->GetPoint( pts[( i+2 )%num_pts], x3.data() );
    vec3_t v2_old = x2 - x0_old;
    vec3_t v3_old = x3 - x0_old;
    
    //top point
    vec3_t S = v2_old.cross( v3_old );
    double V_old = tetraVol( x0_old, S, x2, x3, true );
    double V_new = tetraVol( x0_new, S, x2, x3, true );
    double prod = V_old * V_new;
    if ( prod < 0 ) {
      return ( true );
    }
  }
  return( false );
}
