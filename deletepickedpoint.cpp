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

void DeletePickedPoint::foobar(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, vector <bool> DeadNode)
{
  vtkIdType src_N_points=src->GetNumberOfPoints();
  vector <vtkIdType> OffSet(src_N_points);
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < src_N_points; ++src_id_node) {
    if(!DeadNode[src_id_node])
    {
      vec3_t x;
      src->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(src, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
    else
    {
      //search closest node
//       getsho
        
    }
  };
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
    vtkIdType N_pts, *src_pts, *dst_pts;
    vtkIdType type_cell = src->GetCellType(id_cell);
    src->GetCellPoints(id_cell, N_pts, src_pts);
    src->GetCellPoints(id_cell, N_pts, dst_pts);
    bool DeadCell=false;
    for(int i=0;i<N_pts;i++)
    {
      if(DeadNode[src_pts[i]]) {DeadCell=true;}
      dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
    }
    if(!DeadCell)
    {
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
    }
  };
}

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
  foobar(grid,grid_tmp,hitlist);
  cout_grid(cout,grid_tmp,true,true,true,true);
  makeCopy(grid_tmp, grid);
};
