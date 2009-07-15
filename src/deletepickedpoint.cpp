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

#include "deletepickedpoint.h"

#include <QObject>
#include <QVector>

#include "guimainwindow.h"

DeletePickedPoint::DeletePickedPoint() : RemovePoints()
{
  EG_TYPENAME;
  //Activate undo/redo
  setQuickSave(true);
}

void DeletePickedPoint::operate()
{
  vtkIdType id_node = GuiMainWindow::pointer()->getPickedPoint();
  cout << "You picked " << id_node << endl;

  char type;
  QVector <vtkIdType> PSP;
  
  UpdatePotentialSnapPoints(true);
  
  QMessageBox msgBox;
  msgBox.setText("Delete point?");
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  switch (msgBox.exec()) {
  case QMessageBox::Yes:
    cout<<"yes was clicked"<<endl;
    DeletePoint(id_node);
    break;
  case QMessageBox::No:
    cout<<"no was clicked"<<endl;
    cout<<"=== Topological neighbours ==="<<endl;
    PSP = getPotentialSnapPoints(id_node);
    cout<<"id_node="<<id_node<<" PSP="<<PSP<<endl;
    
    cout<<"=== NODE TYPE ==="<<endl;
    type = getNodeType(id_node);
    cout<<"id_node="<<id_node<<" is of type="<<(int)type<<"="<<VertexType2Str(type)<<endl;
    
    break;
  default:
     // should never be reached
    break;
  }
  
};

bool DeletePickedPoint::DeletePoint( vtkIdType id_node )
{
  int N1 = grid->GetNumberOfPoints();
  
  QVector<vtkIdType> selected_cells;
  getSurfaceCells( m_bcs, selected_cells, grid );
  QVector <vtkIdType> selected_nodes;
  getNodesFromCells( selected_cells, selected_nodes, grid );
  
  setAllSurfaceCells();
  l2l_t n2n = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t nodes = getPartNodes();
  
  UpdatePotentialSnapPoints(true);
  
  EG_VTKDCN( vtkCharArray, node_type, grid, "node_type" );
  EG_VTKDCC( vtkIntArray, cell_code, grid, "cell_code" );
  EG_VTKDCN( vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired" );
  
  // global values
  QSet <vtkIdType> all_deadcells;
  QSet <vtkIdType> all_mutatedcells;
  QSet <vtkIdType> all_mutilatedcells;
  int num_newpoints = 0;
  int num_newcells = 0;
  
  QVector <bool> marked_nodes(nodes.size());
  marked_nodes.fill(false);
  
  QVector <vtkIdType> deadnode_vector;
  QVector <vtkIdType> snappoint_vector;
  
  if ( node_type->GetValue( id_node ) != VTK_FIXED_VERTEX ) {
    if ( !marked_nodes[id_node] && removePointCriteria( id_node ) ) {
      // local values
      QSet <vtkIdType> l_deadcells;
      QSet <vtkIdType> l_mutatedcells;
      QSet <vtkIdType> l_mutilatedcells;
      int l_num_newpoints = 0;
      int l_num_newcells = 0;
      vtkIdType snap_point = FindSnapPoint( id_node, l_deadcells, l_mutatedcells, l_mutilatedcells, l_num_newpoints, l_num_newcells );
      if( snap_point >= 0 ) {
        // add deadnode/snappoint pair
        deadnode_vector.push_back( id_node );
        snappoint_vector.push_back( snap_point );
        // update global values
        num_newpoints += l_num_newpoints;
        num_newcells += l_num_newcells;
        all_deadcells.unite( l_deadcells ); //all_deadcells unite! Kill the living! :D
        all_mutatedcells.unite( l_mutatedcells );
        all_mutilatedcells.unite( l_mutilatedcells );
        // mark neighbour nodes
        foreach( int i_node_neighbour, n2n[_nodes[id_node]] ) marked_nodes[nodes[i_node_neighbour]] = true;
      }
    }
  }
  
  //delete
  DeleteSetOfPoints( deadnode_vector, snappoint_vector, all_deadcells, all_mutatedcells, all_mutilatedcells, num_newpoints, num_newcells );
  
  int N2 = grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
  
  return( m_NumRemoved == 1 );
}
