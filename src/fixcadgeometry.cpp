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
#include "fixcadgeometry.h"
#include "surfacemesher.h"
#include "vertexmeshdensity.h"
#include "guimainwindow.h"

#include <vtkMath.h>

#include "geometrytools.h"
using namespace GeometryTools;

#include <QtDebug>
#include <iostream>
using namespace std;

fixCadGeometry::fixCadGeometry()
{
  EG_TYPENAME;
  m_PerformGeometricTests = true;
  m_UseProjectionForSmoothing = true;
  m_UseNormalCorrectionForSmoothing = true;
  m_FeatureAngle = GeometryTools::deg2rad(200);//this angle is also used by swaptriangles!!!
}

void fixCadGeometry::operate()
{
  qDebug()<<"==>fixing CAD geometry...";
//   prepare();
  setAllCells();
  readSettings();
//   readVMD();
  
  //prepare BCmap
  QSet <int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  QMap <int,int> BCmap;
  foreach(int bc, bcs) BCmap[bc]=1;
  
  //set density infinite
  VertexMeshDensity VMD;
  VMD.density = 9000;
  VMD.BCmap = BCmap;
  cout<<"VMD="<<VMD<<endl;
  qWarning()<<"VMD.BCmap="<<VMD.BCmap;
  m_VMDvector.push_back(VMD);
  cout<<"m_VMDvector="<<m_VMDvector<<endl;
  
  //update node info
//   customUpdateNodeInfo(true);
  updateNodeInfo(true);
  
  //call surface mesher
  setGrid(grid);
  setBoundaryCodes(bcs);
  setVertexMeshDensityVector(m_VMDvector);
//   setDesiredLength();
//   mesher();

  // finalize
  createIndices(grid);
  customUpdateNodeInfo(false);
  setDesiredLength();
}

void fixCadGeometry::mesher()
{
  setDesiredLength();
  if (m_BoundaryCodes.size() == 0) {
    return;
  }
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, grid, "node_meshdensity_desired");
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    characteristic_length_desired->SetValue(id_node, 1e-6);
  }
  int num_inserted = 0;
  int num_deleted = 0;
  int iter = 0;
  bool done = false;
  while (!done) {
    ++iter;
    cout << "surface mesher iteration " << iter << ":" << endl;
    setDesiredLength();
    cout << "  inserted nodes : " << num_inserted << endl;
    customUpdateNodeInfo();
    setDesiredLength();
    int num_deleted = deleteNodes();
    cout << "  deleted nodes  : " << num_deleted << endl;
    swap();
    setDesiredLength();
    int N_crit = grid->GetNumberOfPoints()/100;
    done = (iter >= m_NumMaxIter) || ((num_inserted - num_deleted < N_crit) && (num_inserted + num_deleted < N_crit));
    cout << "  total nodes    : " << grid->GetNumberOfPoints() << endl;
    cout << "  total cells    : " << grid->GetNumberOfCells() << endl;
  }
  createIndices(grid);
  customUpdateNodeInfo(false);
  setDesiredLength();
}

void fixCadGeometry::setDesiredLength(double L)
{
  setAllSurfaceCells();
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2n   = getPartN2N();
  l2l_t  c2c   = getPartC2C();
  
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired,   grid, "node_meshdensity_desired");
  EG_VTKDCN(vtkIntArray,    characteristic_length_specified, grid, "node_specified_density");
  
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    vtkIdType id_node = nodes[i_nodes];
    characteristic_length_specified->SetValue(id_node, 0);
    characteristic_length_desired->SetValue(id_node, L);
  }
}

void fixCadGeometry::customUpdateNodeInfo(bool update_type)
{
  setAllCells();
  l2g_t nodes = getPartNodes();
  foreach (vtkIdType id_node, nodes) {
    if(update_type) {
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      node_type->SetValue(id_node, custom_getNodeType(id_node));
    }
    EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid, "node_meshdensity_current");//what we have
    node_meshdensity_current->SetValue(id_node, CurrentVertexAvgDist(id_node));
    
    EG_VTKDCN(vtkIntArray, node_specified_density, grid, "node_specified_density");//density index from table
    VertexMeshDensity nodeVMD = getVMD(id_node);
    int idx = m_VMDvector.indexOf(nodeVMD);
    node_specified_density->SetValue(id_node, idx);
  }
  writeGrid(grid, "info");
}

char fixCadGeometry::custom_getNodeType( vtkIdType id_node, bool allow_feature_edge_vertices, bool fix_unselected )
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();
  
  //initialize default value
  char type = VTK_FIXED_VERTEX;
  
  //loop through edges around id_node
  
  QVector <vtkIdType> edges;
  
  double CosEdgeAngle = cos(m_EdgeAngle);
  
  foreach( int i_node2, n2n[_nodes[id_node]] ) {
    vtkIdType id_node2 = nodes[i_node2];
    //-----------------------
    //determine edge type
    char edge = custom_getEdgeType( id_node2, id_node, allow_feature_edge_vertices, fix_unselected );
    
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

char fixCadGeometry::custom_getEdgeType( vtkIdType a_node1, vtkIdType a_node2, bool allow_feature_edge_vertices, bool fix_unselected )
{
  double CosFeatureAngle = cos(deg2rad(179));
  
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
    if(m_BoundaryCodes.isEmpty()) {
      EG_BUG;
    }
    if(fix_unselected) {
      if( !m_BoundaryCodes.contains(cell_code_0) || !m_BoundaryCodes.contains(cell_code_1) ) {
        edge = VTK_FIXED_VERTEX;// does not make sense, but should make the points of the edge fixed
      }
    }
  }
  
  return( edge );
}
