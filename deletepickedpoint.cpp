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

void DeletePickedPoint::operate()
{
/*  vtkPointPicker *PointPicker = GuiMainWindow::pointer()->getPointPicker();
  GuiMainWindow::pointer()->getInteractor()->SetPicker(PointPicker);*/
  vtkIdType nodeId = GuiMainWindow::pointer()->getPickedPoint();
  cout<<"You picked "<<nodeId<<endl;
  
  int N_points=grid->GetNumberOfPoints();
  int N_cells=grid->GetNumberOfCells();
  vector <bool> hitlist(N_points);
  vector <vtkIdType> offset(N_points);
  hitlist[nodeId]=true;
  int N_newpoints=-1;
  int N_newcells=0;
  map <vtkIdType,bool> marked;
  foreach(vtkIdType C, n2c[nodeId])
  {
    if(!marked[C]){
      N_newcells-=1;
      marked[C]=true;
    }
  }
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  
  EG_VTKSP(vtkUnstructuredGrid, grid_tmp);
  allocateGrid(grid_tmp,N_cells+N_newcells,N_points+N_newpoints);
  makeCopyNoAllocFiltered(grid,grid_tmp,hitlist);
  cout_grid(cout,grid_tmp,true,true,true,true);
  makeCopy(grid_tmp, grid);
};
