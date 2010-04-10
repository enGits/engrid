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
  if( id_node<0  || GuiMainWindow::pointer()->getPickedObject()!=1 ) {
    QApplication::restoreOverrideCursor();
    EG_ERR_RETURN("Error: No node picked.");
  }
  
  char type;
  QVector <vtkIdType> PSP;
  
/*  QSet <int> bcs;
  GuiMainWindow::pointer()->getAllBoundaryCodes(bcs);
  qWarning()<<"bcs="<<bcs;*/
  
  GuiMainWindow::pointer()->getAllBoundaryCodes(this->m_BoundaryCodes);//IMPORTANT: to make sure only unselected nodes become fixed (redundant with previous line, but more readable)
  qWarning()<<"m_BoundaryCodes="<<m_BoundaryCodes;
  
//   this->setBoundaryCodes(m_BoundaryCodes);
//   qWarning()<<"m_BoundaryCodes="<<m_BoundaryCodes;
  
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

bool DeletePickedPoint::DeletePoint(vtkIdType id_node)
{
  int N1 = m_Grid->GetNumberOfPoints();
  
  QVector<vtkIdType> selected_cells;
  getSurfaceCells(m_BoundaryCodes, selected_cells, m_Grid);
  QVector<vtkIdType> selected_nodes;
  getNodesFromCells(selected_cells, selected_nodes, m_Grid);
  
  setAllSurfaceCells();
  l2l_t  n2n   = getPartN2N();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  nodes = getPartNodes();
  
  UpdatePotentialSnapPoints(false);
  
  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type" );
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code" );
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired" );
  
  // global values
  QVector <vtkIdType> all_deadcells;
  QVector <vtkIdType> all_mutatedcells;
  int num_newpoints = 0;
  int num_newcells = 0;
  
  QVector <bool> marked_nodes(nodes.size(), false);
  
  QVector <vtkIdType> deadnode_vector;
  QVector <vtkIdType> snappoint_vector;
  
  if (node_type->GetValue(id_node) != VTK_FIXED_VERTEX) {
      // local values
      QVector <vtkIdType> dead_cells;
      QVector <vtkIdType> mutated_cells;
      int l_num_newpoints = 0;
      int l_num_newcells = 0;
    
      setDebugLevel(11);
      vtkIdType snap_point = FindSnapPoint(id_node, dead_cells, mutated_cells, l_num_newpoints, l_num_newcells, marked_nodes);
      qDebug()<<"snap_point="<<snap_point;
      setDebugLevel(0);
    
      if(snap_point >= 0) {
        // add deadnode/snappoint pair
        deadnode_vector.push_back(id_node);
        snappoint_vector.push_back(snap_point);
        // update global values
        num_newpoints += l_num_newpoints;
        num_newcells  += l_num_newcells;
        all_deadcells += dead_cells;
        all_mutatedcells += mutated_cells;
        // mark neighbour nodes
        foreach(int i_node_neighbour, n2n[_nodes[id_node]]) {
          marked_nodes[nodes[i_node_neighbour]] = true;
        }
      }
    }
  
  //delete
  DeleteSetOfPoints(deadnode_vector, snappoint_vector, all_deadcells, all_mutatedcells, num_newpoints, num_newcells);
  
  int N2 = m_Grid->GetNumberOfPoints();
  m_NumRemoved = N1 - N2;
  
  return( m_NumRemoved == 1 );
}
