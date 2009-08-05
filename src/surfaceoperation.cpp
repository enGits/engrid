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
#include "surfaceoperation.h"

#include "guimainwindow.h"

#include <vtkCharArray.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkPolygon.h>

#include "geometrytools.h"
using namespace GeometryTools;

SurfaceOperation::SurfaceOperation()
    : Operation()
{
  //default values for determining node types and for smoothing operations
  m_Convergence = 0;
  m_NumberOfIterations = 20;
  m_RelaxationFactor = 0.01;
  m_AllowFeatureEdgeVertices = 1;//0 by default in VTK, but we need 1 to avoid the "potatoe effect" ^^
  m_FeatureAngle = GeometryTools::deg2rad(45);
  setFeatureAngle(GeometryTools::deg2rad(45));
  getSet("surface meshing", "edge angle to determine fixed vertices", 180, m_EdgeAngle);
  m_EdgeAngle = GeometryTools::deg2rad(m_EdgeAngle);
  setEdgeAngle(m_EdgeAngle);
  m_BoundarySmoothing = 1;
}

void SurfaceOperation::operate()
{

}

ostream& operator<<( ostream &out, stencil_t S )
{
  out << "S.id_cell1=" << S.id_cell1 << " ";
  out << "S.id_cell2=" << S.id_cell2 << " ";
  out << "S.sameBC=" << S.sameBC << " ";
  out << "S.twocells=" << S.twocells << " ";
  out << "S.neighbour_type=" << S.neighbour_type << " ";
  out << "[";
  for ( int i = 0; i < 4; i++ ) {
    out << S.p[i];
    if ( i != 3 ) out << ",";
  }
  out << "]";
  return( out );
}

stencil_t SurfaceOperation::getStencil( vtkIdType id_cell1, int j1 )
{
  l2g_t  cells = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2l_t  c2c   = getPartC2C();

  if ( grid->GetCellType( id_cell1 ) != VTK_TRIANGLE ) {
    cout << "CELL IS NOT A TRIANGLE" << endl;
    EG_BUG;
  }

  //return variable
  stencil_t S;

  //default values:
  S.sameBC = false;
  S.twocells = false;
  S.neighbour_type = -1;

  //initialize first cell
  S.id_cell1 = id_cell1;
  vtkIdType N1, *pts1;
  grid->GetCellPoints( S.id_cell1, N1, pts1 );
  //place points 0,1,3
  if ( j1 == 0 ) { S.p[0] = pts1[2]; S.p[1] = pts1[0]; S.p[3] = pts1[1]; }
  else if ( j1 == 1 ) { S.p[0] = pts1[0]; S.p[1] = pts1[1]; S.p[3] = pts1[2]; }
  else if ( j1 == 2 ) { S.p[0] = pts1[1]; S.p[1] = pts1[2]; S.p[3] = pts1[0]; };

  //initialize second cell
  S.id_cell2 = -1;
  S.p[2] = -1;

  //twocells
  if ( c2c[_cells[id_cell1]][j1] != -1 ) { //if neighbour cell

    //twocells
    S.twocells = true;
    S.id_cell2 = cells[c2c[_cells[id_cell1]][j1]];

    //sameBC
    EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
    if ( cell_code->GetValue( S.id_cell1 ) == cell_code->GetValue( S.id_cell2 ) ) S.sameBC = true;

    //neighbour_type
    S.neighbour_type = grid->GetCellType( S.id_cell2 );
    if ( S.neighbour_type == VTK_TRIANGLE ) {//if neighbour cell is a triangle
      vtkIdType N2, *pts2;
      grid->GetCellPoints( S.id_cell2, N2, pts2 );

      //place point 2
      bool p2 = false;
      if ( c2c[_cells[S.id_cell2]][0] != -1 ) {
        if ( cells[c2c[_cells[S.id_cell2]][0]] == S.id_cell1 ) {
          S.p[2] = pts2[2];
          p2 = true;
        }
      }
      if ( c2c[_cells[S.id_cell2]][1] != -1 ) {
        if ( cells[c2c[_cells[S.id_cell2]][1]] == S.id_cell1 ) {
          S.p[2] = pts2[0];
          p2 = true;
        }
      }
      if ( c2c[_cells[S.id_cell2]][2] != -1 ) {
        if ( cells[c2c[_cells[S.id_cell2]][2]] == S.id_cell1 ) {
          S.p[2] = pts2[1];
          p2 = true;
        }
      }

      if ( !p2 ) {//failed to place point 2, appears when cell1 is linked to cell2, but cell2 not to cell1
        cout << "S.id_cell1=" << S.id_cell1 << endl;
        cout << "S.id_cell2=" << S.id_cell2 << endl;
        GuiMainWindow::pointer()->saveAs( GuiMainWindow::pointer()->getFilePath() + "abort.egc", false );
        EG_BUG;
      }
    }
  }//end of if neighbour cell
  return S;
}

int SurfaceOperation::UpdateCurrentMeshDensity()
{
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity START===" << endl;
  }
  QVector<vtkIdType> cells;
  getAllSurfaceCells( cells, grid );
  EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );
  setGrid( grid );
  setCells( cells );
  if ( DebugLevel > 5 ) {
    cout << "cells.size()=" << cells.size() << endl;
  }
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current" );
  l2g_t nodes = getPartNodes();
  foreach( vtkIdType node, nodes ) {
    node_meshdensity_current->SetValue( node, CurrentMeshDensity( node ) );
  }
  if ( DebugLevel > 0 ) {
    cout << "===UpdateMeshDensity END===" << endl;
  }
  return( 0 ); ///@@@ what for???
}

int SurfaceOperation::UpdatePotentialSnapPoints( bool update_node_types, bool allow_feature_edge_vertices )
{
  setAllSurfaceCells();

  l2g_t nodes  = getPartNodes();
  l2g_t cells  = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2l_t c2c    = getPartC2C();

  m_PotentialSnapPoints.resize( nodes.size() );

  //initialize default values
  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  foreach( vtkIdType id_node, nodes ) {
    if ( update_node_types ) node_type->SetValue( id_node, VTK_SIMPLE_VERTEX );
    m_PotentialSnapPoints[id_node].clear();
  }

  //cout<<"===pre-processing==="<<endl;
  int num_edges = 0;
  //We loop through edges
  foreach( vtkIdType id_cell, cells ) {
    vtkIdType *pts, Npts;
    grid->GetCellPoints( id_cell, Npts, pts );
    for ( int i = 0; i < Npts; i++ ) {

      int i_neighbour_cell = c2c[_cells[id_cell]][i];
      if ( i_neighbour_cell >= 0 && cells[i_neighbour_cell] < id_cell ) continue;//already visited edge
      num_edges++;

      vtkIdType id_node1 = pts[i];
      vtkIdType id_node2 = pts[( i+1 )%Npts];

      //-----------------------
      //determine edge type
      char edge = getEdgeType( id_node2, id_node1, allow_feature_edge_vertices );
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
        grid->GetPoint( m_PotentialSnapPoints[id_node][0], x1 );
        grid->GetPoint( id_node, x2 );
        grid->GetPoint( m_PotentialSnapPoints[id_node][1], x3 );
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

char SurfaceOperation::getNodeType( vtkIdType id_node, bool allow_feature_edge_vertices )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  //initialize default value
  char type = VTK_SIMPLE_VERTEX;

  //loop through edges around id_node

  QVector <vtkIdType> edges;

  double CosEdgeAngle = cos(this->m_EdgeAngle) ;

  foreach( int i_node2, n2n[_nodes[id_node]] ) {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = getEdgeType( id_node2, id_node, allow_feature_edge_vertices );

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
      grid->GetPoint( edges[0], x1 );
      grid->GetPoint( id_node, x2 );
      grid->GetPoint( edges[1], x3 );
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

  if ( !allow_feature_edge_vertices && type == VTK_FEATURE_EDGE_VERTEX ) EG_BUG;
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

char SurfaceOperation::getEdgeType( vtkIdType a_node1, vtkIdType a_node2, bool allow_feature_edge_vertices, bool fix_unselected )
{
  double CosFeatureAngle = cos(this->m_FeatureAngle);

  //compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int numNei = getEdgeCells( a_node1, a_node2, neighbour_cells ) - 1;

  //set default value
  char edge = VTK_SIMPLE_VERTEX;

  if ( numNei == 0 ) {
    edge = VTK_BOUNDARY_EDGE_VERTEX;
  }
  else if ( numNei >= 2 ) {
    qWarning() << "FATAL ERROR: edge belongs to more than 2 cells! This is not supported yet.";
    EG_BUG;
    edge = VTK_FEATURE_EDGE_VERTEX;
  }
  else if ( numNei == 1 ) {
    //check angle between cell1 and cell2 against FeatureAngle
    if ( allow_feature_edge_vertices && this->m_AllowFeatureEdgeVertices && CosAngle( grid, neighbour_cells[0], neighbour_cells[1] ) <= CosFeatureAngle ) {
      edge = VTK_FEATURE_EDGE_VERTEX;
    }
    //check the boundary codes
    EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
    int cell_code_0 = cell_code->GetValue( neighbour_cells[0] );
    int cell_code_1 = cell_code->GetValue( neighbour_cells[1] );
    if ( cell_code_0 !=  cell_code_1 ) {
      edge = VTK_BOUNDARY_EDGE_VERTEX;
    }
//     qWarning()<<"m_BoundaryCodes="<<m_BoundaryCodes;
    if(m_BoundaryCodes.isEmpty()) abort();
    if(fix_unselected) {
      if( !m_BoundaryCodes.contains(cell_code_0) || !m_BoundaryCodes.contains(cell_code_1) ) {
        edge = VTK_FIXED_VERTEX;// does not make sense, but should make the points of the edge fixed
      }
    }
  }

  return( edge );
}

QSet <int> SurfaceOperation::getBCset( vtkIdType id_node )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();

  EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
  QSet <int> bc;
  foreach( int i_cell, n2c[_nodes[id_node]] ) {
    vtkIdType id_cell = cells[i_cell];
    bc.insert( cell_code->GetValue( id_cell ) );
  }
  return( bc );
}

VertexMeshDensity SurfaceOperation::getVMD( vtkIdType id_node )
{
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();

  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );

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
double SurfaceOperation::CurrentVertexAvgDist( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_dist = 0;
  double avg_dist = 0;
  int N = n2n[_nodes[id_node]].size();
  vec3_t C;
  grid->GetPoint( id_node, C.data() );
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    grid->GetPoint( id_node_neighbour, M.data() );
    total_dist += ( M - C ).abs();
  }
  avg_dist = total_dist / ( double )N;
  return( avg_dist );
}

double SurfaceOperation::CurrentMeshDensity( vtkIdType id_node )
{
  return 1.0 / CurrentVertexAvgDist( id_node );
}

double SurfaceOperation::DesiredVertexAvgDist( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_dist = 0;
  double avg_dist = 0;
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );
  int N = n2n[_nodes[id_node]].size();
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    total_dist += 1. / node_meshdensity_desired->GetValue( id_node_neighbour );
  }
  avg_dist = total_dist / ( double )N;
  return( avg_dist );
}

double SurfaceOperation::DesiredMeshDensity( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  double total_density = 0;
  double avg_density = 0;
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );
  int N = n2n[_nodes[id_node]].size();
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    total_density += node_meshdensity_desired->GetValue( id_node_neighbour );
  }
  avg_density = total_density / ( double )N;
  return( avg_density );
}

///////////////////////////////////////////

//---------------------------------------------------
//Utility functions used in Roland's formulas

///@@@ TODO: change meshdensity fields to edgelength fields since this is what is mostly used?

/// desired edge length for id_node
double SurfaceOperation::desiredEdgeLength( vtkIdType id_node )
{
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );
  return( 1.0 / node_meshdensity_desired->GetValue( id_node ) );
}

//other functions
///perimeter
double SurfaceOperation::perimeter( vtkIdType id_cell )
{
  double ret = 0;
  vtkIdType num_pts, *pts;
  grid->GetCellPoints( id_cell, num_pts, pts );
  for ( int i = 0; i < num_pts; i++ ) {
    vec3_t A, B;
    grid->GetPoints()->GetPoint( pts[i], A.data() );
    grid->GetPoints()->GetPoint( pts[( i+1 )%num_pts], B.data() );
    ret += ( B - A ).abs();
  }
  return( ret );
}

/// mean desired edge length for id_cell
double SurfaceOperation::meanDesiredEdgeLength( vtkIdType id_cell )
{
  vtkIdType num_pts, *pts;
  grid->GetCellPoints( id_cell, num_pts, pts );
  int total = 0;
  for ( int i = 0; i < num_pts; i++ ) {
    total += desiredEdgeLength( pts[i] );
  }
  return total / ( double )num_pts;
}

///@@@ TODO: Should be renamed to be more explicit if possible

/// perimeter / sum of the desired edge lengths
double SurfaceOperation::Q_L( vtkIdType id_cell )
{
  double denom_sum = 0;
  vtkIdType num_pts, *pts;
  grid->GetCellPoints( id_cell, num_pts, pts );
  for ( int i = 0; i < num_pts; i++ ) {
    denom_sum += desiredEdgeLength( pts[i] );
  }
  return( perimeter( id_cell ) / denom_sum );
}

/// sum(2*edgelength,edges(id_node))/sum(desired edgelengths of each edgepoint,edges(id_node))
double SurfaceOperation::Q_L1( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes = getPartNodes();

  double num_sum = 0;
  double denom_sum = 0;
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    num_sum += 2 * distance( grid, id_node_neighbour, id_node );
    denom_sum += desiredEdgeLength( id_node ) + desiredEdgeLength( id_node_neighbour );
  }
  return( num_sum / denom_sum );
}

/// minimum of sum(2*edgelength)/sum(desired edgelengths of each edgepoint) for each edge of id_node
double SurfaceOperation::Q_L2( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes = getPartNodes();

  QVector <double> V;
  double num, denom;
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    num = 2 * distance( grid, id_node_neighbour, id_node );
    denom = desiredEdgeLength( id_node ) + desiredEdgeLength( id_node_neighbour );
    V.push_back( num / denom );
  }
  qSort( V.begin(), V.end() );
  return( V[0] );
}

/// Value to minimize for mesh smoothing. w allows putting more weight on the form or the area of triangles.
double SurfaceOperation::T_min( int w )
{
  l2g_t cells = getPartCells();
  double T = 0;
  foreach( vtkIdType id_cell, cells ) {
    T += areaOfCircumscribedCircle( grid, id_cell ) / pow( cellVA( grid, id_cell ), w ) * pow( meanDesiredEdgeLength( id_cell ), 2 * ( w - 1 ) );
  }
  return( T );
}

//---------------------------------------------------

vtkIdType SurfaceOperation::getClosestNode( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes = getPartNodes();

  vec3_t C;
  grid->GetPoint( id_node, C.data() );
  vtkIdType id_minlen = -1;
  double minlen = -1;
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    grid->GetPoint( id_node_neighbour, M.data() );
    double len = ( M - C ).abs();
    if ( minlen < 0 or len < minlen ) {
      minlen = len;
      id_minlen = id_node_neighbour;
    }
  }
  return( id_minlen );
}

vtkIdType SurfaceOperation::getFarthestNode( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes = getPartNodes();

  vec3_t C;
  grid->GetPoint( id_node, C.data() );
  vtkIdType id_maxlen = -1;
  double maxlen = -1;
  foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) {
    vtkIdType id_node_neighbour = nodes[i_node_neighbour];
    vec3_t M;
    grid->GetPoint( id_node_neighbour, M.data() );
    double len = ( M - C ).abs();
    if ( maxlen < 0 or len > maxlen ) {
      maxlen = len;
      id_maxlen = id_node_neighbour;
    }
  }
  return( id_maxlen );
}

QVector <vtkIdType> SurfaceOperation::getPotentialSnapPoints( vtkIdType id_node )
{
  if ((id_node < 0) || (id_node >= m_PotentialSnapPoints.size())) {
    EG_BUG;
  }
  return m_PotentialSnapPoints[id_node];
}
