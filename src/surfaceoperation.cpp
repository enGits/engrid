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
  Convergence = 0;
  NumberOfIterations = 20;
  RelaxationFactor = 0.01;
  FeatureEdgeSmoothing = 1;//0 by default in VTK, but we need 1 to avoid the "potatoe effect" ^^
  FeatureAngle = 45;
  EdgeAngle = 15;
  BoundarySmoothing = 1;
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

int SurfaceOperation::UpdateNodeType()
{
  l2g_t nodes  = getPartNodes();
  l2g_t cells  = getPartCells();
  g2l_t _cells = getPartLocalCells();
  l2l_t c2c    = getPartC2C();

  //cout<<"=== UpdateNodeType START ==="<<endl;
  //prepare
  setAllSurfaceCells();

  m_PotentialSnapPoints.resize( nodes.size() );

  //initialize default values
  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  foreach( vtkIdType id_node, nodes ) {
    node_type->SetValue( id_node, VTK_SIMPLE_VERTEX );
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
      char edge = getEdgeType( id_node2, id_node1 );
      //-----------------------
      //determine node type pre-processing (count nb of complex edges if the node is complex, otherwise, just count the nb of edges)
      if ( edge && node_type->GetValue( id_node1 ) == VTK_SIMPLE_VERTEX ) {
        m_PotentialSnapPoints[id_node1].clear();
        m_PotentialSnapPoints[id_node1].push_back( id_node2 );
        node_type->SetValue( id_node1, edge );
      }
      else if (( edge && node_type->GetValue( id_node1 ) == VTK_BOUNDARY_EDGE_VERTEX ) ||
               ( edge && node_type->GetValue( id_node1 ) == VTK_FEATURE_EDGE_VERTEX ) ||
               ( !edge && node_type->GetValue( id_node1 ) == VTK_SIMPLE_VERTEX ) ) {
        m_PotentialSnapPoints[id_node1].push_back( id_node2 );
        if ( node_type->GetValue( id_node1 ) && edge == VTK_BOUNDARY_EDGE_VERTEX ) {
          node_type->SetValue( id_node1, VTK_BOUNDARY_EDGE_VERTEX );//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }

      if ( edge && node_type->GetValue( id_node2 ) == VTK_SIMPLE_VERTEX ) {
        m_PotentialSnapPoints[id_node2].clear();
        m_PotentialSnapPoints[id_node2].push_back( id_node1 );
        node_type->SetValue( id_node2, edge );
      }
      else if (( edge && node_type->GetValue( id_node2 ) == VTK_BOUNDARY_EDGE_VERTEX ) ||
               ( edge && node_type->GetValue( id_node2 ) == VTK_FEATURE_EDGE_VERTEX ) ||
               ( !edge && node_type->GetValue( id_node2 ) == VTK_SIMPLE_VERTEX ) ) {
        m_PotentialSnapPoints[id_node2].push_back( id_node1 );
        if ( node_type->GetValue( id_node2 ) && edge == VTK_BOUNDARY_EDGE_VERTEX ) {
          node_type->SetValue( id_node2, VTK_BOUNDARY_EDGE_VERTEX );//VTK_BOUNDARY_EDGE_VERTEX has priority over VTK_FEATURE_EDGE_VERTEX
        }
      }
    }
  }

  //cout<<"num_edges="<<num_edges<<endl;

  //-----------------------
  //determine node type post-processing
  double CosEdgeAngle = cos(( double ) vtkMath::RadiansFromDegrees( this->EdgeAngle ) );
  //cout<<"===post-processing==="<<endl;
  //This time, we loop through nodes
  foreach( vtkIdType id_node, nodes ) {
    if ( node_type->GetValue( id_node ) == VTK_FEATURE_EDGE_VERTEX || node_type->GetValue( id_node ) == VTK_BOUNDARY_EDGE_VERTEX ) { //see how many edges; if two, what the angle is

      if ( !this->BoundarySmoothing &&
           node_type->GetValue( id_node ) == VTK_BOUNDARY_EDGE_VERTEX ) {
        node_type->SetValue( id_node, VTK_FIXED_VERTEX );
      }
      else if ( m_PotentialSnapPoints[id_node].size() != 2 ) {
        node_type->SetValue( id_node, VTK_FIXED_VERTEX );
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
          node_type->SetValue( id_node, VTK_FIXED_VERTEX );
        }
      }//if along edge
    }//if edge vertex
  }
  //cout<<"m_PotentialSnapPoints.size()="<<m_PotentialSnapPoints.size()<<endl;
  //cout<<"=== UpdateNodeType END ==="<<endl;
  return( 0 );
}

///@@@  TODO: Optimize
char SurfaceOperation::getNodeType( vtkIdType id_node )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  //initialize default value
  char type = VTK_SIMPLE_VERTEX;

  //loop through edges around id_node

  QVector <vtkIdType> edges;

  double CosEdgeAngle = cos(( double ) vtkMath::RadiansFromDegrees( this->EdgeAngle ) );

  foreach( int i_node2, n2n[_nodes[id_node]] ) {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = getEdgeType( id_node2, id_node );

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

    if ( !this->BoundarySmoothing &&
         type == VTK_BOUNDARY_EDGE_VERTEX ) {
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

char SurfaceOperation::getEdgeType( vtkIdType a_node1, vtkIdType a_node2 )
{
  double CosFeatureAngle = cos(( double ) vtkMath::RadiansFromDegrees( this->FeatureAngle ) );

  //compute number of cells around edge [a_node,p2] and put them into neighbour_cells
  QVector <vtkIdType> neighbour_cells;
  int numNei = getEdgeCells( a_node1, a_node2, neighbour_cells ) - 1;

  //set default value
  char edge = VTK_SIMPLE_VERTEX;

  if ( numNei == 0 ) {
    edge = VTK_BOUNDARY_EDGE_VERTEX;
  }
  else if ( numNei >= 2 ) {
    edge = VTK_FEATURE_EDGE_VERTEX;
  }
  else if ( numNei == 1 ) {
    //check angle between cell1 and cell2 against FeatureAngle
    if ( this->FeatureEdgeSmoothing && CosAngle( grid, neighbour_cells[0], neighbour_cells[1] ) <= CosFeatureAngle ) {
      edge = VTK_FEATURE_EDGE_VERTEX;
    }
    //check the boundary codes
    EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
    if ( cell_code->GetValue( neighbour_cells[0] ) !=  cell_code->GetValue( neighbour_cells[1] ) ) {
      edge = VTK_BOUNDARY_EDGE_VERTEX;
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

///@@@ TODO: Check that operations using n2n,n2c,c2c are still correct
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
  double num_sum = 0;
  double denom_sum = 0;
  foreach( vtkIdType j, n2n[id_node] ) {
    num_sum += 2 * distance( grid, j, id_node );
    denom_sum += desiredEdgeLength( id_node ) + desiredEdgeLength( j );
  }
  return( num_sum / denom_sum );
}

/// minimum of sum(2*edgelength)/sum(desired edgelengths of each edgepoint) for each edge of id_node
double SurfaceOperation::Q_L2( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  QVector <double> V;
  double num, denom;
  foreach( vtkIdType j, n2n[id_node] ) {
    num = 2 * distance( grid, j, id_node );
    denom = desiredEdgeLength( id_node ) + desiredEdgeLength( j );
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
  vec3_t C;
  grid->GetPoint( id_node, C.data() );
  vtkIdType id_minlen = -1;
  double minlen = -1;
  foreach( vtkIdType neighbour, n2n[id_node] ) {
    vec3_t M;
    grid->GetPoint( neighbour, M.data() );
    double len = ( M - C ).abs();
    if ( minlen < 0 or len < minlen ) {
      minlen = len;
      id_minlen = neighbour;
    }
  }
  return( id_minlen );
}

vtkIdType SurfaceOperation::getFarthestNode( vtkIdType id_node )
{
  l2l_t n2n = getPartN2N();
  vec3_t C;
  grid->GetPoint( id_node, C.data() );
  vtkIdType id_maxlen = -1;
  double maxlen = -1;
  foreach( vtkIdType neighbour, n2n[id_node] ) {
    vec3_t M;
    grid->GetPoint( neighbour, M.data() );
    double len = ( M - C ).abs();
    if ( maxlen < 0 or len > maxlen ) {
      maxlen = len;
      id_maxlen = neighbour;
    }
  }
  return( id_maxlen );
}

int SurfaceOperation::NumberOfCommonPoints( vtkIdType id_node1, vtkIdType id_node2, bool& IsTetra )
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
    // TEST 0: id_node1, id_node2 and intersection* must form a cell
    QVector <vtkIdType> EdgeCells_1i;
    QVector <vtkIdType> EdgeCells_2i;
    QVector <vtkIdType> inter;
    int N;
    
    // intersection1
    N = getEdgeCells( id_node1, intersection1, EdgeCells_1i );
    if(N!=2) EG_BUG;
    N = getEdgeCells( id_node2, intersection1, EdgeCells_2i );
    if(N!=2) EG_BUG;
    qcontIntersection(EdgeCells_1i, EdgeCells_2i, inter);
    if(inter.size()<=0) EG_BUG;
    
    // intersection2
    N = getEdgeCells( id_node1, intersection2, EdgeCells_1i );
    if(N!=2) EG_BUG;
    N = getEdgeCells( id_node2, intersection2, EdgeCells_2i );
    if(N!=2) EG_BUG;
    qcontIntersection(EdgeCells_1i, EdgeCells_2i, inter);
    if(inter.size()<=0) EG_BUG;
    
    // TEST 1
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

QVector <vtkIdType> SurfaceOperation::getPotentialSnapPoints( vtkIdType id_node )
{
  if ( id_node < 0 || id_node >= m_PotentialSnapPoints.size() ) EG_BUG;
  return m_PotentialSnapPoints[id_node];
}

// DEFINITIONS:
// Normal cell: nothing has changed
// Dead cell: the cell does not exist anymore
// Mutated cell: the cell's form has changed
// Mutilated cell: the cell has less points than before

///@@@  TODO: Organize cases and make sure all are considered if possible.
vtkIdType SurfaceOperation::FindSnapPoint( vtkIdType DeadNode, QSet <vtkIdType> & DeadCells, QSet <vtkIdType> & MutatedCells, QSet <vtkIdType> & MutilatedCells, int& num_newpoints, int& num_newcells )
{
  // preparations
  setAllSurfaceCells();
  l2l_t n2c = getPartN2C();

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
    foreach( vtkIdType id_cell, n2c[DeadNode] ) { //loop through potentially dead cells
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
        if ( pts[i] != DeadNode && pts[i] != PSP &&  n2c[pts[i]].size() <= 1 ) {
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

        // TEST 4
        if ( Old_N*New_N<0 || New_N*New_N < Old_N*Old_N*1. / 100. ) { //area + inversion check
          if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
          IsValidSnapPoint = false;
        }
        
        ///@@@ TODO: finish this
        // TEST 5: flipped cell test from old laplace smoother
//         if(FlippedCells( DeadNode, vec3_t P )) {
//           qWarning()<<"EPIC FAIL!!!!!!!!!!!!!!!!!! You have flipped cells!";
//           EG_BUG;
//           if ( DebugLevel > 10 ) cout << "Sorry, but you are not allowed to move point " << DeadNode << " to point " << PSP << "." << endl;
//           IsValidSnapPoint = false;
//         }

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

bool SurfaceOperation::DeletePoint( vtkIdType DeadNode, int& num_newpoints, int& num_newcells )
{
  QSet <vtkIdType> DeadNodes;
  DeadNodes.insert( DeadNode );
  bool ret = DeleteSetOfPoints( DeadNodes, num_newpoints, num_newcells );
  return( ret );
}
//End of DeletePoint

bool SurfaceOperation::DeleteSetOfPoints( QSet <vtkIdType> DeadNodes, int& num_newpoints, int& num_newcells )
{
  int initial_num_points = grid->GetNumberOfPoints();
  
  CheckSurfaceIntegrity check_surface_integrity;
  check_surface_integrity();
  if(!check_surface_integrity.isWaterTight()) {
    qWarning()<<"FATAL ERROR: NOT WATERTIGHT!";
    GuiMainWindow::pointer()->saveAs( GuiMainWindow::pointer()->getFilePath() + "abort.egc", false );
    EG_BUG;
  }
  
  QVector<vtkIdType> deadnode_vector = Set2Vector( DeadNodes, false );

//   UpdateNodeType();

  //src grid info
  int num_points = grid->GetNumberOfPoints();
  int num_cells = grid->GetNumberOfCells();

  QSet <vtkIdType> all_deadcells;
  QSet <vtkIdType> all_mutatedcells;
  QSet <vtkIdType> all_mutilatedcells;
  QVector <vtkIdType> SnapPoint( deadnode_vector.size() );

  //counter init
  num_newpoints = 0;
  num_newcells = 0;

  for ( int i = 0; i < deadnode_vector.size(); i++ ) {
    if ( deadnode_vector[i] < 0 || deadnode_vector[i] >= num_points ) {
      cout << "Warning: Point out of range: deadnode_vector[i]=" << deadnode_vector[i] << " num_points=" << num_points << endl;
      EG_BUG;
      return ( false );
    }

    //local values
    int l_num_newpoints = 0;
    int l_num_newcells = 0;
    QSet <vtkIdType> l_DeadCells;
    QSet <vtkIdType> l_MutatedCells;
    QSet <vtkIdType> l_MutilatedCells;

    SnapPoint[i] = FindSnapPoint( deadnode_vector[i], l_DeadCells, l_MutatedCells, l_MutilatedCells, l_num_newpoints, l_num_newcells );
    //global values
    num_newpoints += l_num_newpoints;
    num_newcells += l_num_newcells;
    all_deadcells.unite( l_DeadCells ); //all_deadcells unite! Kill the living! :D
    all_mutatedcells.unite( l_MutatedCells );
    all_mutilatedcells.unite( l_MutilatedCells );

    if ( DebugLevel > 0 ) {
      cout << "===>deadnode_vector[i]=" << deadnode_vector[i] << " moving to SNAPPOINT=" << SnapPoint[i] << " DebugLevel=" << DebugLevel << endl;
    }
    if ( SnapPoint[i] < 0 ) {
      cout << "ERROR: no possible SnapPoint found." << endl;
      EG_BUG;
      return( false );
    }

  }

  if(num_newcells!=2*num_newpoints) {
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
//       dst_pts = new vtkIdType[dst_num_pts];
      
      if ( all_mutatedcells.contains( id_cell ) ) { //mutated cell
        int num_deadnode = 0;
        for ( int i = 0; i < src_num_pts; i++ ) {
          int DeadIndex = deadnode_vector.indexOf( src_pts[i] );
          if ( DeadIndex != -1 ) {
            dst_pts[i] = SnapPoint[DeadIndex] - OffSet[SnapPoint[DeadIndex]]; // dead node
            num_deadnode++;
          }
          else {
            dst_pts[i] = src_pts[i] - OffSet[src_pts[i]]; // not a dead node
          }
        }
        if(num_deadnode!=1) {
          qWarning()<<"FATAL ERROR: Mutated cell has more than one dead node!";
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
            qWarning()<<"FATAL ERROR: Normal cell contains a dead node!";
            EG_BUG;
          }
          dst_pts[i] = src_pts[i] - OffSet[src_pts[i]];
        }
      }
      // copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell( type_cell, dst_num_pts, dst_pts );
      copyCellData( grid, id_cell, dst, id_new_cell );
//       delete dst_pts;
    }
  }

  CheckSurfaceIntegrity check_surface_integrity_tmp;
  check_surface_integrity_tmp.setGrid(dst);
  check_surface_integrity_tmp();
  if(!check_surface_integrity_tmp.isWaterTight()) {
    qWarning()<<"FATAL ERROR: NOT WATERTIGHT!";
    GuiMainWindow::pointer()->saveAs( GuiMainWindow::pointer()->getFilePath() + "pre_abort.egc", false );
    makeCopy( dst, grid );
    GuiMainWindow::pointer()->saveAs( GuiMainWindow::pointer()->getFilePath() + "abort.egc", false );
    int final_num_points = grid->GetNumberOfPoints();
    if ( initial_num_points - final_num_points != DeadNodes.size() ) {
      EG_BUG;
    }
    EG_BUG;
  }
  
  makeCopy( dst, grid );
  int final_num_points = grid->GetNumberOfPoints();
  
  if ( -num_newpoints != DeadNodes.size() ) {
    EG_BUG;
  }

  if ( initial_num_points - final_num_points != DeadNodes.size() ) {
    EG_BUG;
  }
  
  check_surface_integrity();
  if(!check_surface_integrity.isWaterTight()) {
    qWarning()<<"FATAL ERROR: NOT WATERTIGHT!";
    GuiMainWindow::pointer()->saveAs( GuiMainWindow::pointer()->getFilePath() + "abort.egc", false );
    EG_BUG;
  }
  
  return( true );
}
//End of DeleteSetOfPoints

bool SurfaceOperation::FlippedCells( vtkIdType id_node, vec3_t P )
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
