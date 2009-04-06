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
#include "guimainwindow.h"
#include "egvtkobject.h"
#include "geometrytools.h"
using namespace GeometryTools;

void DeletePickedPoint::operate()
{
  vtkIdType nodeId = GuiMainWindow::pointer()->getPickedPoint();
  cout<<"You picked "<<nodeId<<endl;

//   DeletePoint(grid,nodeId);
  vtkIdType Boss, Peon1, Peon2;
  int BC=0;
  Boss=nodeId;
  getNeighbours(Boss,Peon1,Peon2,BC);
  cout<<"Boss="<<Boss<<" Peon1="<<Peon1<<" Peon2="<<Peon2<<" BC="<<BC<<endl;
};
