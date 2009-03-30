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

void DeletePickedPoint::foobar(vtkUnstructuredGrid *src, vtkUnstructuredGrid *dst, vtkIdType DeadNode)
{
  //Find closest point to DeadNode
  vtkIdType SnapPoint = getClosestNode(DeadNode,src);//DeadNode moves to SnapPoint
  cout<<"SnapPoint="<<SnapPoint<<endl;
  
  //src grid info
  int N_points=src->GetNumberOfPoints();
  int N_cells=src->GetNumberOfCells();
  
  //vector used to redefine the new point IDs
  vector <vtkIdType> OffSet(N_points);
  
  //count number of points and cells to remove + allocate + put mutated cells in a set
  int N_newpoints=-1;
  int N_newcells=0;
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  foreach(vtkIdType C, n2c[DeadNode])//loop through potentially dead cells
  {
    //get points around cell
    vtkIdType N_pts, *pts;
    src->GetCellPoints(C, N_pts, pts);
    
    bool ContainsSnapPoint=false;
    for(int i=0;i<N_pts;i++)
    {
      cout<<"pts["<<i<<"]="<<pts[i]<<" and SnapPoint="<<SnapPoint<<endl;
      if(pts[i]==SnapPoint) {ContainsSnapPoint=true;break;}
//       newcell[i]= ( (pts[i]==src_id_node) ? SnapPoint : pts[i] );
    }
    if(ContainsSnapPoint)
    {
      if(N_pts<=3)//dead cell
      {
        DeadCells.insert(C);
        N_newcells-=1;
        cout<<"cell "<<C<<" has been pwned!"<<endl;
      }
      else//mutilated cell
      {
        MutilatedCells.insert(C);
        cout<<"cell "<<C<<" has lost a limb!"<<endl;
      }
    }
    else
    {
      //mutated cell
      MutatedCells.insert(C);
      cout<<"cell "<<C<<" has been infected!"<<endl;
    }
  }
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  
  //copy undead points
  vtkIdType dst_id_node=0;
  for (vtkIdType src_id_node = 0; src_id_node < N_points; src_id_node++) {//loop through src points
    if(src_id_node!=DeadNode)//if the node isn't dead, copy it
    {
      vec3_t x;
      src->GetPoints()->GetPoint(src_id_node, x.data());
      dst->GetPoints()->SetPoint(dst_id_node, x.data());
      copyNodeData(src, src_id_node, dst, dst_id_node);
      OffSet[src_id_node]=src_id_node-dst_id_node;
      dst_id_node++;
    }
  };
  
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType N_pts, *src_pts, *dst_pts;
      vtkIdType type_cell = src->GetCellType(id_cell);
      src->GetCellPoints(id_cell, N_pts, src_pts);
      src->GetCellPoints(id_cell, N_pts, dst_pts);
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        for(int i=0;i<N_pts;i++)
        {
          if(src_pts[i]==DeadNode) dst_pts[i]=SnapPoint-OffSet[SnapPoint];
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      if(MutilatedCells.contains(id_cell))//mutilated cell
      {
        if(type_cell==VTK_QUAD) type_cell==VTK_TRIANGLE;
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        int j=0;
        for(int i=0;i<N_pts;i++)
        {
          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead
          //do nothing in case of DeadNode
        }
      }
      else//normal cell
      {
        for(int i=0;i<N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
    }
  };
}

void DeletePickedPoint::operate()
{
  vtkIdType nodeId = GuiMainWindow::pointer()->getPickedPoint();
  cout<<"You picked "<<nodeId<<endl;
  
  EG_VTKSP(vtkUnstructuredGrid, grid_tmp);
  foobar(grid,grid_tmp,nodeId);
  cout_grid(cout,grid_tmp,true,true,true,true);
  makeCopy(grid_tmp, grid);
};
