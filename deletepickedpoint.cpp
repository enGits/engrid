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

// DEFINITIONS:
// Normal cell: nothing has changed
// Dead cell: the cell does not exist anymore
// Mutated cell: the cell's form has changed
// Mutilated cell: the cell has less points than before

int DeletePickedPoint::NumberOfCommonPoints(vtkIdType node1, vtkIdType node2, bool& IsTetra)
{
//   QVector< QSet< int > > 	n2n
  QSet <int> node1_neighbours=n2n[node1];
  QSet <int> node2_neighbours=n2n[node2];
  QSet <int> intersection=node1_neighbours.intersect(node2_neighbours);
  int N=intersection.size();
  IsTetra=false;
  if(N==2)
  {
    QSet<int>::const_iterator p1=intersection.begin();
    QSet<int>::const_iterator p2=p1+1;
    cout<<"*p1="<<*p1<<endl;
    cout<<"*p2="<<*p2<<endl;
    vtkIdType intersection1=_nodes[*p1];
    vtkIdType intersection2=_nodes[*p2];
    if(n2n[intersection1].contains(intersection2))//if there's an edge between intersection1 and intersection2
    {
      //check if (node1,intersection1,intersection2) and (node2,intersection1,intersection2) are defined as cells!
  //     QVector< QSet< int > > 	n2c
      QSet< int > S1=n2c[intersection1];
      QSet< int > S2=n2c[intersection2];
      QSet< int > Si=S1.intersect(S2);
      int counter=0;
      foreach(vtkIdType C,Si){
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(C, N_pts, pts);
        for(int i=0;i<N_pts;i++)
        {
          if(pts[i]==node1 || pts[i]==node2) counter++;
        }
      }
      cout<<"counter="<<counter<<endl;
      if(counter>=2) IsTetra=true;
    }
  }
  return(N);
}

vtkIdType DeletePickedPoint::FindSnapPoint(vtkUnstructuredGrid *src, vtkIdType DeadNode)
{

}

bool DeletePickedPoint::DeletePoint(vtkUnstructuredGrid *src, vtkIdType DeadNode)
{
  EG_VTKSP(vtkUnstructuredGrid, dst);
  
    //src grid info
  int N_points=src->GetNumberOfPoints();
  int N_cells=src->GetNumberOfCells();
  int N_newpoints=-1;
  int N_newcells=0;
  
  QSet <vtkIdType> DeadCells;
  QSet <vtkIdType> MutatedCells;
  QSet <vtkIdType> MutilatedCells;
  
  vtkIdType SnapPoint=-1;
  //Find closest point to DeadNode
//   vtkIdType SnapPoint = getClosestNode(DeadNode,src);//DeadNode moves to SnapPoint
  
  foreach(vtkIdType PSP, n2n[DeadNode])
  {
    bool IsValidSnapPoint=true;
    
    cout<<"====>PSP="<<PSP<<endl;
    bool IsTetra=true;
    if(NumberOfCommonPoints(DeadNode,PSP,IsTetra)>2)//common point check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    if(IsTetra)//tetra check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    
    //count number of points and cells to remove + analyse cell transformations
    N_newpoints=-1;
    N_newcells=0;
    DeadCells.clear();
    MutatedCells.clear();
    MutilatedCells.clear();
    foreach(vtkIdType C, n2c[DeadNode])//loop through potentially dead cells
    {
      //get points around cell
      vtkIdType N_pts, *pts;
      src->GetCellPoints(C, N_pts, pts);
      
      bool ContainsSnapPoint=false;
      for(int i=0;i<N_pts;i++)
      {
        cout<<"pts["<<i<<"]="<<pts[i]<<" and PSP="<<PSP<<endl;
        if(pts[i]==PSP) {ContainsSnapPoint=true;break;}
//         if(n2c[pts[i]]<=1) invincible=true;
      }
      if(ContainsSnapPoint)
      {
        if(N_pts==3)//dead cell
        {
          //TODO: Check that empty lines aren't left behind when a cell is killed
//           if(invincible)
          DeadCells.insert(C);
          N_newcells-=1;
          cout<<"cell "<<C<<" has been pwned!"<<endl;
        }
  /*      else if(N_pts==4)//mutilated cell
        {
          MutilatedCells.insert(C);
          cout<<"cell "<<C<<" has lost a limb!"<<endl;
        }*/
        else
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
      }
      else
      {
        vtkIdType src_N_pts, *src_pts;
        src->GetCellPoints(C, src_N_pts, src_pts);
        
        if(src_N_pts!=3)
        {
          cout<<"RED ALERT: Xenomorph detected!"<<endl;
          EG_BUG;
        }
        
        vtkIdType OldTriangle[3];
        vtkIdType NewTriangle[3];
        
        for(int i=0;i<src_N_pts;i++)
        {
          OldTriangle[i]=src_pts[i];
          NewTriangle[i]=( (src_pts[i]==DeadNode) ? PSP : src_pts[i] );
        }
        vec3_t Old_N= triNormal(src, OldTriangle[0], OldTriangle[1], OldTriangle[2]);
        vec3_t New_N= triNormal(src, NewTriangle[0], NewTriangle[1], NewTriangle[2]);
        double OldArea=Old_N.abs();
        double NewArea=New_N.abs();
        double scal=Old_N*New_N;
        double cross=(Old_N.cross(New_N)).abs();//double-cross on Nar Shadaa B-)
        
        cout<<"OldArea="<<OldArea<<endl;
        cout<<"NewArea="<<NewArea<<endl;
        cout<<"scal="<<scal<<endl;
        cout<<"cross="<<cross<<endl;
        
        if(Old_N*New_N<Old_N*Old_N*1./100.)//area + inversion check
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
  /*      if(NewArea<OldArea*1./100.)
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }
        
        if(abs(cross)>10e-4)
        {
          cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
          IsValidSnapPoint=false;
        }*/
        
        //mutated cell
        MutatedCells.insert(C);
        cout<<"cell "<<C<<" has been infected!"<<endl;
      }
    }
    
    if(N_cells+N_newcells<=0)//survivor check
    {
      cout<<"Sorry, but you are not allowed to move point "<<DeadNode<<" to point "<<PSP<<"."<<endl;
      IsValidSnapPoint=false;
    }
    if(IsValidSnapPoint) {SnapPoint=PSP; break;}
  }//end of loop through potential SnapPoints
  
  cout<<"===>SNAPPOINT="<<SnapPoint<<endl;
  if(SnapPoint<0) {cout<<"Sorry no possible SnapPoint found."<<endl; return(false);}
  
  //allocate
  cout<<"N_points="<<N_points<<endl;
  cout<<"N_cells="<<N_cells<<endl;
  cout<<"N_newpoints="<<N_newpoints<<endl;
  cout<<"N_newcells="<<N_newcells<<endl;
  allocateGrid(dst,N_cells+N_newcells,N_points+N_newpoints);
  
  //vector used to redefine the new point IDs
  QVector <vtkIdType> OffSet(N_points);
  
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
  
  cout<<"DeadCells="<<DeadCells<<endl;
  cout<<"MutatedCells="<<MutatedCells<<endl;
  cout<<"MutilatedCells="<<MutilatedCells<<endl;
  
  //Copy undead cells
  for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {//loop through src cells
    if(!DeadCells.contains(id_cell))//if the cell isn't dead
    {
      vtkIdType src_N_pts, *src_pts;
      vtkIdType dst_N_pts, *dst_pts;
      src->GetCellPoints(id_cell, src_N_pts, src_pts);
      
      vtkIdType type_cell = src->GetCellType(id_cell);
      cout<<"-->id_cell="<<id_cell<<endl;
      for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
//       src->GetCellPoints(id_cell, dst_N_pts, dst_pts);
      dst_N_pts=src_N_pts;
      dst_pts=new vtkIdType[dst_N_pts];
      if(MutatedCells.contains(id_cell))//mutated cell
      {
        cout<<"processing mutated cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==DeadNode) {
            cout<<"SnapPoint="<<SnapPoint<<endl;
            cout<<"OffSet[SnapPoint]="<<OffSet[SnapPoint]<<endl;
            cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
            dst_pts[i]=SnapPoint-OffSet[SnapPoint];
          }
          else dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
        cout<<"--->dst_pts:"<<endl;
        for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
        
      }
      else if(MutilatedCells.contains(id_cell))//mutilated cell
      {
        cout<<"processing mutilated cell "<<id_cell<<endl;
        
        if(type_cell==VTK_QUAD) {
          type_cell=VTK_TRIANGLE;
          dst_N_pts-=1;
        }
        else {cout<<"FATAL ERROR: Unknown mutilated cell detected! It is not a quad! Potential xenomorph infestation!"<<endl;EG_BUG;}
        //merge points
        int j=0;
        for(int i=0;i<src_N_pts;i++)
        {
          if(src_pts[i]==SnapPoint) { dst_pts[j]=SnapPoint-OffSet[SnapPoint];j++; }//SnapPoint
          else if(src_pts[i]!=DeadNode) { dst_pts[j]=src_pts[i]-OffSet[src_pts[i]];j++; }//pre-snap/dead + post-snap/dead
          //do nothing in case of DeadNode
        }
      }
      else//normal cell
      {
        cout<<"processing normal cell "<<id_cell<<endl;
        for(int i=0;i<src_N_pts;i++)
        {
          dst_pts[i]=src_pts[i]-OffSet[src_pts[i]];
        }
      }
      //copy the cell
      vtkIdType id_new_cell = dst->InsertNextCell(type_cell, dst_N_pts, dst_pts);
      copyCellData(src, id_cell, dst, id_new_cell);
      cout<<"===Copying cell "<<id_cell<<" to "<<id_new_cell<<"==="<<endl;
      cout<<"src_pts:"<<endl;
      for(int i=0;i<src_N_pts;i++) cout<<"src_pts["<<i<<"]="<<src_pts[i]<<endl;
      cout<<"dst_pts:"<<endl;
      for(int i=0;i<dst_N_pts;i++) cout<<"dst_pts["<<i<<"]="<<dst_pts[i]<<endl;
      cout<<"OffSet="<<OffSet<<endl;
      cout<<"===Copying cell end==="<<endl;
      delete dst_pts;
    }
  };
  cout_grid(cout,dst,true,true,true,true);
  makeCopy(dst, src);
  return(true);
}

void DeletePickedPoint::operate()
{
  vtkIdType nodeId = GuiMainWindow::pointer()->getPickedPoint();
  cout<<"You picked "<<nodeId<<endl;

  DeletePoint(grid,nodeId);
};
