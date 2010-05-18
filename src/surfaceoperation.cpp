//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2010 enGits GmbH                                     +
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
#include "surfaceoperation.h"

#include "guimainwindow.h"

#include <vtkCharArray.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>

#include "geometrytools.h"
using namespace GeometryTools;

SurfaceOperation::SurfaceOperation() : Operation()
{
  //default values for determining node types and for smoothing operations
  getSet("surface meshing", "edge angle to determine fixed vertices", 180, m_EdgeAngle);
  getSet("surface meshing", "feature angle", 180, m_FeatureAngle);
  m_FeatureAngle = GeometryTools::deg2rad(m_FeatureAngle);
  m_EdgeAngle = GeometryTools::deg2rad(m_EdgeAngle);
  setEdgeAngle(m_EdgeAngle);
  m_BoundarySmoothing = 1;
}

void SurfaceOperation::operate()
{

}

ostream& operator<<(ostream &out, stencil_t S)
{
  out << "S.id_cell = " << S.id_cell << " ";
  out << "S.id_node = " << S.id_node << " ";
  out << "S.sameBC = " << S.sameBC << " ";
  out << "S.type = " << S.type_cell << " ";
  out << "S.p1 = " << S.p1 << " ";
  out << "S.p2 = " << S.p2 << " ";
  return(out);
}

stencil_t SurfaceOperation::getStencil(vtkIdType id_cell1, int j1)
{
  stencil_t S;
  {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell1, N_pts, pts);
    S.p1 = pts[j1];
    S.p2 = pts[0];
    if (j1 < N_pts - 1) {
      S.p2 = pts[j1 + 1];
    }
  }
  QSet<vtkIdType> cells_p1;
  for (int i = 0; i < m_Part.n2cGSize(S.p1); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(S.p1, i);
    if (id_cell != id_cell1) {
      cells_p1.insert(id_cell);
    }
  }
  QSet<vtkIdType> cells_p2;
  for (int i = 0; i < m_Part.n2cGSize(S.p2); ++i) {
    vtkIdType id_cell = m_Part.n2cGG(S.p2, i);
    if (id_cell != id_cell1) {
      cells_p2.insert(id_cell);
    }
  }
  QSet<vtkIdType> cells = cells_p1.intersect(cells_p2);
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  S.sameBC = true;
  S.id_cell.resize(1);
  S.id_cell[0] = id_cell1;
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      S.id_cell.push_back(id_cell);
      if (cell_code->GetValue(id_cell) != cell_code->GetValue(id_cell1)) {
        S.sameBC = false;
      }
    }
  }
  S.id_node.resize(S.id_cell.size());
  S.type_cell.resize(S.id_cell.size());
  for (int i = 0; i < S.id_cell.size(); ++i) {
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(S.id_cell[i], N_pts, pts);
    S.type_cell[i] = m_Grid->GetCellType(S.id_cell[i]);
    for (int j = 0; j < N_pts; ++j) {
      if (pts[j] != S.p1 && pts[j] != S.p2) {
        S.id_node[i] = pts[j];
        break;
      }
    }
  }
  return S;
}

int SurfaceOperation::UpdateCurrentMeshDensity()
{
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity START===" << endl;
  }
  QVector<vtkIdType> cells;
  getAllSurfaceCells( cells, m_Grid );
  EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );
  EG_VTKDCN( vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired" );
  setGrid( m_Grid );
  setCells( cells );
  if ( DebugLevel > 5 ) {
    cout << "cells.size()=" << cells.size() << endl;
  }
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_current, m_Grid, "node_meshdensity_current" );
  l2g_t nodes = getPartNodes();
  foreach( vtkIdType node, nodes ) {
    node_meshdensity_current->SetValue( node, CurrentMeshDensity( node ) );
  }
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity END===" << endl;
  }
  return( 0 ); ///\todo what for???
}

int SurfaceOperation::UpdatePotentialSnapPoints( bool update_node_types, bool fix_unselected)
{
  setAllSurfaceCells();

  l2g_t nodes  = getPartNodes();
  l2g_t cells  = getPartCells();

  m_PotentialSnapPoints.resize(m_Grid->GetNumberOfPoints());

  //initialize default values
  EG_VTKDCN( vtkCharArray, node_type, m_Grid, "node_type" );
  foreach( vtkIdType id_node, nodes ) {
    if ( update_node_types ) node_type->SetValue( id_node, VTK_SIMPLE_VERTEX );
    m_PotentialSnapPoints[id_node].clear();
  }

  //cout<<"===pre-processing==="<<endl;
  int num_edges = 0;
  //We loop through edges
  foreach( vtkIdType id_cell, cells ) {
    vtkIdType *pts, Npts;
    m_Grid->GetCellPoints( id_cell, Npts, pts );
    for ( int i = 0; i < Npts; i++ ) {

      QVector <vtkIdType> EdgeCells;
      int Ncells = getEdgeCells( pts[i], pts[(i+1)%Npts], EdgeCells );
      
      bool already_visited_edge = false;
      for(int i_cell=0;i_cell<Ncells;i_cell++) {
        if(EdgeCells[i_cell]<id_cell) {
          already_visited_edge = true;
          break;
        }
      }
      
      if ( already_visited_edge ) continue;//already visited edge
      num_edges++;

      vtkIdType id_node1 = pts[i];
      vtkIdType id_node2 = pts[( i+1 )%Npts];

      //-----------------------
      //determine edge type
      char edge = getEdgeType( id_node2, id_node1, fix_unselected );
      //-----------------------
      //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
      if ( edge && node_type->GetValue( id_node1 ) == VTK_SIMPLE_VERTEX ) {
        m_PotentialSnapPoints[id_node1].clear();
        m_PotentialSnapPoints[id_node1].push_back( id_node2 );
        if ( update_node_types ) node_type->SetValue( id_node1, edge );
      }
      else if (( edge && node_type->GetValue( id_node1 ) == VTK_BOUNDARY_EDGE_VERTEX ) ||
               ( edge && node_type->GetValue( id_node1 ) == VTK_FEATURE_EDGE_VERTEX ) ||
               ( !edge && node_type->GetValue( id_node1 ) == VTK_SIMPLE_VERTEX ) ) {
        m_PotentialSnapPoints[id_node1].push_back( id_node2 );
        if ( node_type->GetValue( id_node1 ) && edge == VTK_BOUNDARY_EDGE_VERTEX ) {
          if ( update_node_types ) node_type->SetValue( id_node1, VTK_BOUNDARY_EDGE_VERTEX );//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }

      if ( edge && node_type->GetValue( id_node2 ) == VTK_SIMPLE_VERTEX ) {
        m_PotentialSnapPoints[id_node2].clear();
        m_PotentialSnapPoints[id_node2].push_back( id_node1 );
        if ( update_node_types ) node_type->SetValue( id_node2, edge );
      }
      else if (( edge && node_type->GetValue( id_node2 ) == VTK_BOUNDARY_EDGE_VERTEX ) ||
               ( edge && node_type->GetValue( id_node2 ) == VTK_FEATURE_EDGE_VERTEX ) ||
               ( !edge && node_type->GetValue( id_node2 ) == VTK_SIMPLE_VERTEX ) ) {
        m_PotentialSnapPoints[id_node2].push_back( id_node1 );
        if ( node_type->GetValue( id_node2 ) && edge == VTK_BOUNDARY_EDGE_VERTEX ) {
          if ( update_node_types ) node_type->SetValue( id_node2, VTK_BOUNDARY_EDGE_VERTEX );//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }
    }
  }

  //cout<<"num_edges="<<num_edges<<endl;

  //-----------------------
  //determine node type post-processing
  double CosEdgeAngle = cos(this->m_EdgeAngle);
  //cout<<"===post-processing==="<<endl;
  //This time, we loop through nodes
  foreach( vtkIdType id_node, nodes ) {
    if ( node_type->GetValue( id_node ) == VTK_FEATURE_EDGE_VERTEX || node_type->GetValue( id_node ) == VTK_BOUNDARY_EDGE_VERTEX ) { //see how many edges; if two, what the angle is

      if ( !this->m_BoundarySmoothing && node_type->GetValue( id_node ) == VTK_BOUNDARY_EDGE_VERTEX ) {
        if ( update_node_types ) node_type->SetValue( id_node, VTK_FIXED_VERTEX );
      } else if ( m_PotentialSnapPoints[id_node].size() != 2 ) {
        if ( update_node_types ) node_type->SetValue( id_node, VTK_FIXED_VERTEX );
      }
      else { //check angle between edges
        double x1[3], x2[3], x3[3], l1[3], l2[3];
        m_Grid->GetPoint( m_PotentialSnapPoints[id_node][0], x1 );
        m_Grid->GetPoint( id_node, x2 );
        m_Grid->GetPoint( m_PotentialSnapPoints[id_node][1], x3 );
        for ( int k = 0; k < 3; k++ ) {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if ( vtkMath::Normalize( l1 ) >= 0.0 &&
             vtkMath::Normalize( l2 ) >= 0.0 &&
             vtkMath::Dot( l1, l2 ) < CosEdgeAngle ) {
          if ( update_node_types ) node_type->SetValue( id_node, VTK_FIXED_VERTEX );
        }
      }//if along edge
    }//if edge vertex
  }
  //cout<<"m_PotentialSnapPoints.size()="<<m_PotentialSnapPoints.size()<<endl;
  //cout<<"=== UpdatePotentialSnapPoints END ==="<<endl;
  return( 0 );
}

char SurfaceOperation::getNodeType( vtkIdType id_node, bool fix_unselected )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  //initialize default value
  char type = VTK_SIMPLE_VERTEX;

  //loop through edges around id_node

  QVector <vtkIdType> edges;

  double CosEdgeAngle = cos(this->m_EdgeAngle);

  foreach( int i_node2, n2n[_nodes[id_node]] ) {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = getEdgeType(id_node2, id_node, fix_unselected);

    //-----------------------
    //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
    if ( edge && type == VTK_SIMPLE_VERTEX ) {
      edges.clear();
      edges.push_back( id_node2 );
      type = edge;
    }
    else if (( edge && type == VTK_BOUNDARY_EDGE_VERTEX ) ||
             ( edge && type == VTK_FEATURE_EDGE_VERTEX ) ||
             ( !edge && type == VTK_SIMPLE_VERTEX ) ) {
      edges.push_back( id_node2 );
      if ( type && edge == VTK_BOUNDARY_EDGE_VERTEX ) {
        type = VTK_BOUNDARY_EDGE_VERTEX;//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
      }
    }
  }
  //-----------------------
  //determine node type post-processing
  if ( type == VTK_FEATURE_EDGE_VERTEX || type == VTK_BOUNDARY_EDGE_VERTEX ) { //see how many edges; if two, what the angle is

    if ( !this->m_BoundarySmoothing && type == VTK_BOUNDARY_EDGE_VERTEX ) {
      type = VTK_FIXED_VERTEX;
    }
    else if ( edges.size() != 2 ) {
      type = VTK_FIXED_VERTEX;
    }
    else { //check angle between edges
      double x1[3], x2[3], x3[3], l1[3], l2[3];
      m_Grid->GetPoint( edges[0], x1 );
      m_Grid->GetPoint( id_node, x2 );
      m_Grid->GetPoint( edges[1], x3 );
      for ( int k = 0; k < 3; k++ ) {
        l1[k] = x2[k] - x1[k];
        l2[k] = x3[k] - x2[k];
      }
      if ( vtkMath::Normalize( l1 ) >= 0.0 &&
           vtkMath::Normalize( l2 ) >= 0.0 &&
           vtkMath::Dot( l1, l2 ) < CosEdgeAngle ) {
             type = VTK_FIXED_VERTEX;
      }
    }//if along edge
  }//if edge vertex

  return( type );
}

int SurfaceOperation::getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QVector <vtkIdType> &EdgeCells )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t cells  = getPartCells();
  l2l_t n2c    = getPartN2C();

  QSet<vtkIdType> S1;
  foreach( int i, n2c[_nodes[id_node1]] ) {
    S1.insert( cells[i] );
  }

  QSet<vtkIdType> S2;
  foreach( int i, n2c[_nodes[id_node2]] ) {
    S2.insert( cells[i] );
  }

  S2.intersect( S1 );
  EdgeCells = Set2Vector( S2, false );
  return EdgeCells.size();
}

int SurfaceOperation::getEdgeCells( vtkIdType id_node1, vtkIdType id_node2, QSet <vtkIdType> &EdgeCells )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t cells  = getPartCells();
  l2l_t n2c    = getPartN2C();

  QSet<vtkIdType> S1;
  foreach( int i, n2c[_nodes[id_node1]] ) {
    S1.insert( cells[i] );
  }

  QSet<vtkIdType> S2;
  foreach( int i, n2c[_nodes[id_node2]] ) {
    S2.insert( cells[i] );
  }

  EdgeCells = S2.intersect( S1 );
  return EdgeCells.size();
}

char SurfaceOperation::getEdgeType(vtkIdType a_node1, vtkIdType a_node2, bool fix_unselected)
{
  double CosFeatureAngle = cos(this->m_FeatureAngle);
  bool feature_edges_disabled = m_FeatureAngle >= M_PI;

  //compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int numNei = getEdgeCells( a_node1, a_node2, neighbour_cells ) - 1;

  //set default value
  char edge = VTK_SIMPLE_VERTEX;

  if ( numNei == 0 ) {
    edge = VTK_BOUNDARY_EDGE_VERTEX;
  }
  else if ( numNei >= 2 ) {
    //qWarning() << "FATAL ERROR: edge belongs to more than 2 cells! This is not supported yet.";
    //EG_BUG;
    //edge = VTK_FEATURE_EDGE_VERTEX;
    edge = VTK_BOUNDARY_EDGE_VERTEX;
  }
  else if ( numNei == 1 ) {
    //check angle between cell1 and cell2 against FeatureAngle
    if (CosAngle(m_Grid, neighbour_cells[0], neighbour_cells[1] ) <= CosFeatureAngle && !feature_edges_disabled) {
      edge = VTK_FEATURE_EDGE_VERTEX;
    }
    //check the boundary codes
    EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );
    int cell_code_0 = cell_code->GetValue( neighbour_cells[0] );
    int cell_code_1 = cell_code->GetValue( neighbour_cells[1] );
    if ( cell_code_0 !=  cell_code_1 ) {
      edge = VTK_BOUNDARY_EDGE_VERTEX;
    }
//     qWarning()<<"m_BoundaryCodes="<<m_BoundaryCodes;
    if(m_BoundaryCodes.isEmpty()) {
      EG_ERR_RETURN("no boundary codes specified");
    }
    if(fix_unselected) {
      if( !m_BoundaryCodes.contains(cell_code_0) || !m_BoundaryCodes.contains(cell_code_1) ) {
        edge = VTK_FIXED_VERTEX;// does not make sense, but should make the points of the edge fixed
      }
    }
  }

  return( edge );
}

VertexMeshDensity SurfaceOperation::getVMD( vtkIdType id_node )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();

  EG_VTKDCN( vtkCharArray, node_type, m_Grid, "node_type" );
  EG_VTKDCC( vtkIntArray, cell_code, m_Grid, "cell_code" );

  VertexMeshDensity VMD;
  VMD.type = node_type->GetValue( id_node );
  VMD.density = 0;
  VMD.CurrentNode = id_node;

  foreach( int i_cell, n2c[_nodes[id_node]] ) {
    vtkIdType id_cell = cells[i_cell];
    VMD.BCmap[cell_code->GetValue( id_cell )] = 2;
  }
  return( VMD );
}

//////////////////////////////////////////////
double SurfaceOperation::currentVertexAvgDist( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_dist = 0;
  double avg_dist = 0;
  int N = n2n[_nodes[id_node]].size();
  vec3_t C;
  m_Grid->GetPoint( id_node, C.data() );
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    m_Grid->GetPoint( id_node_neighbour, M.data() );
    total_dist += ( M - C ).abs();
  }
  avg_dist = total_dist / ( double )N;
  return( avg_dist );
}

double SurfaceOperation::CurrentMeshDensity( vtkIdType id_node )
{
  return 1.0 / currentVertexAvgDist( id_node );
}

///\todo change meshdensity fields to edgelength fields since this is what is mostly used?

/// desired edge length for id_node
double SurfaceOperation::desiredEdgeLength( vtkIdType id_node )
{
  EG_VTKDCN( vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired" );
  return( 1.0 / characteristic_length_desired->GetValue( id_node ) );
}

/// mean desired edge length for id_cell
double SurfaceOperation::meanDesiredEdgeLength( vtkIdType id_cell )
{
  vtkIdType num_pts, *pts;
  m_Grid->GetCellPoints( id_cell, num_pts, pts );
  int total = 0;
  for ( int i = 0; i < num_pts; i++ ) {
    total += desiredEdgeLength( pts[i] );
  }
  return total / ( double )num_pts;
}

QVector <vtkIdType> SurfaceOperation::getPotentialSnapPoints( vtkIdType id_node )
{
  if ((id_node < 0) || (id_node >= m_PotentialSnapPoints.size())) {
    // UpdatePotentialSnapPoints should probably be called before using this function.
    EG_BUG;
  }
  return m_PotentialSnapPoints[id_node];
}

bool SurfaceOperation::isCell(vtkIdType id_node1, vtkIdType id_node2, vtkIdType id_node3)
{
  QVector <vtkIdType> EdgeCells_12;
  QVector <vtkIdType> EdgeCells_13;
  QVector <vtkIdType> inter;
  
  getEdgeCells( id_node1, id_node2, EdgeCells_12 );
  getEdgeCells( id_node1, id_node3, EdgeCells_13 );
  qcontIntersection( EdgeCells_12, EdgeCells_13, inter );
  if(inter.size()>1) {
    qWarning()<<"(id_node1, id_node2, id_node3)="<<"("<<id_node1<<", "<<id_node2<<", "<<id_node3<<")";
    qWarning()<<"EdgeCells_12="<<EdgeCells_12;
    qWarning()<<"EdgeCells_13="<<EdgeCells_13;
    qWarning()<<"inter="<<inter;
    writeGrid(m_Grid, "abort");
    EG_BUG;// multiple cells in the same place
  }
  if(DebugLevel>100) {
    qDebug()<<"(id_node1, id_node2, id_node3)="<<"("<<id_node1<<", "<<id_node2<<", "<<id_node3<<")";
    qDebug()<<"EdgeCells_12="<<EdgeCells_12;
    qDebug()<<"EdgeCells_13="<<EdgeCells_13;
    qDebug()<<"inter="<<inter;
  }
  return(inter.size()>0);
}

vec3_t SurfaceOperation::transform(vtkIdType id_node, vec3_t x)
{
  vec3_t x0;
  m_Grid->GetPoint(id_node, x0.data());
  return x;
}
