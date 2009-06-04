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
#include "egvtkobject.h"
#include "geometrytools.h"
using namespace GeometryTools;

DeletePickedPoint::DeletePickedPoint() : SurfaceOperation()
{
  //Activate undo/redo
  setQuickSave(true);
}

void DeletePickedPoint::operate()
{
  vtkIdType nodeId = GuiMainWindow::pointer()->getPickedPoint();
  cout<<"You picked "<<nodeId<<endl;

  int N_newpoints;
  int N_newcells;
  
  vtkIdType id_node;
  char type;
  QVector <vtkIdType> PSP;
  
  UpdateNodeType();
  
  QMessageBox msgBox;
  msgBox.setText("Delete point?");
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  switch (msgBox.exec()) {
  case QMessageBox::Yes:
    cout<<"yes was clicked"<<endl;
    DeletePoint(nodeId,N_newpoints,N_newcells);
    break;
  case QMessageBox::No:
    cout<<"no was clicked"<<endl;
    id_node=nodeId;
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
