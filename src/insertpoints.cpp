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
#include "insertpoints.h"

#include "guimainwindow.h"

#include <vtkCharArray.h>

#include <QTime>

InsertPoints::InsertPoints()
: SurfaceOperation()
{
  setQuickSave(true);
}

void InsertPoints::operate()
{
  if(insert_FP) insert_FP_all();
  if(insert_EP) insert_EP_all();
}

bool InsertPoints::insert_fieldpoint(vtkIdType D)
{
  double Fred1=1.0/sqrt(3);
  double Qmin=1.1;//1.189;
  double total=0;
  for(int i=0;i<3;i++)
  {
    vtkIdType cell=DN(i,D);
    if(cell!=-1) total += Q_L(cell);
  }
  return ( Q_L(D)>1.0/Fred1 && total>3*Qmin );
}

bool InsertPoints::insert_edgepoint(vtkIdType j,vtkIdType K)// node1 K, node2 j
{
  bool result = L_k(j,K)>0.5*(G_k(j)+G_k(K));
  return ( result );
}

bool InsertPoints::SplitSide(vtkIdType id_cell,int side)
{
  vtkIdType N_pts,*pts;
  grid->GetCellPoints(id_cell,N_pts,pts);
  return( insert_edgepoint(pts[side],pts[(side+1)%N_pts]) );
}

int InsertPoints::insert_FP_all()
{
  cout<<"===insert_FP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  setAllSurfaceCells();
  UpdateNodeType();
  
  QVector <vtkIdType> l_SelectedCells;
  getSurfaceCells(m_bcs, l_SelectedCells, grid);
  
  QVector <bool> l_marked_cells(l_SelectedCells.size());
  
  int l_N_newpoints=0;
  int l_N_newcells=0;
  
  //counter
  for(int i_cell=0;i_cell<l_SelectedCells.size();i_cell++)
  {
    vtkIdType id_cell = l_SelectedCells[i_cell];
    if( insert_fieldpoint(id_cell) )
    {
      l_marked_cells[i_cell] = true;
      l_N_newcells += 2;
      l_N_newpoints += 1;
    }
  }
  
  //initialize grid_tmp
  int l_N_points = grid->GetNumberOfPoints();
  int l_N_cells = grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,l_N_cells+l_N_newcells,l_N_points+l_N_newpoints);
  makeCopyNoAlloc(grid, grid_tmp);
  
  //initialize new node counter
  vtkIdType l_newNodeId = l_N_points;

  //actor
  for(int i_cell=0;i_cell<l_SelectedCells.size();i_cell++)
  {
    vtkIdType id_cell = l_SelectedCells[i_cell];
    if( l_marked_cells[i_cell] )
    {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t C(0,0,0);
      for(int i=0;i<N_pts;i++)
      {
        vec3_t corner;
        grid->GetPoints()->GetPoint(pts[i], corner.data());
        C+=corner;
      }
      C=(1/(double)N_pts)*C;
      
      C=project(C);
      grid_tmp->GetPoints()->SetPoint(l_newNodeId,C.data());
      copyNodeData(grid_tmp,pts[0],grid_tmp,l_newNodeId);
      EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");
      node_type->SetValue(l_newNodeId, VTK_SIMPLE_VERTEX);
      
      for(int i=0;i<N_pts;i++)
      {
        vtkIdType pts_triangle[3];
        pts_triangle[0]=pts[i];
        pts_triangle[1]=pts[(i+1)%N_pts];
        pts_triangle[2]=l_newNodeId;
        if(i==0)
        {
          grid_tmp->ReplaceCell(id_cell , 3, pts_triangle);
        }
        else
        {
          vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle);
          copyCellData(grid_tmp,id_cell,grid_tmp,newCellId);
        }
      }
      l_newNodeId++;
    }
  }
  
  //update grid
  makeCopy(grid_tmp,grid);
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_FP_all END==="<<endl;
  return(0);
}

int InsertPoints::insert_EP_all()
{
  l2g_t  cells = getPartCells();
  g2l_t _cells = getPartLocalCells();

  cout<<"===insert_EP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  setAllSurfaceCells();
  UpdateNodeType();
  
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  int l_N_newpoints=0;
  int l_N_newcells=0;
  
  QVector <int> l_marked_cells(cells.size());
  QVector <stencil_t> l_StencilVector(cells.size());
  
  //counter
  for(int i=0;i<cells.size();i++) {
    vtkIdType id_cell=cells[i];
    if(m_bcs.contains(cell_code->GetValue(id_cell)) && grid->GetCellType(id_cell)==VTK_TRIANGLE) {//if selected and triangle cell
      for(int j = 0; j < 3; ++j) {
        stencil_t S = getStencil(id_cell,j);
        if(S.twocells && S.neighbour_type==VTK_TRIANGLE) {
          if( l_marked_cells[i]==0 && l_marked_cells[_cells[S.id_cell2]]==0 ) {
            if( SplitSide(id_cell,j)) {
              l_StencilVector[i]=S;
              l_marked_cells[i]=1;
              l_marked_cells[_cells[S.id_cell2]]=2;
              l_N_newpoints++;
              l_N_newcells+=2;
            }//end of if SplitSide
          }//end of if unmarked
        }//end of if 2 triangles
        else if(!S.twocells) {
          if( l_marked_cells[i]==0 ) {
            if( SplitSide(id_cell,j)) {
                l_StencilVector[i]=S;
                l_marked_cells[i]=1;
                l_N_newpoints++;
                l_N_newcells+=1;
            }//end of if SplitSide
          }//end of if unmarked
        }//end of if 1 triangle
      }//end of loop through sides
    }//end of if selected and triangle cell
  }//end of counter loop
  
  //initialize grid_tmp
  int l_N_points=grid->GetNumberOfPoints();
  int l_N_cells=grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,l_N_cells+l_N_newcells,l_N_points+l_N_newpoints);
  makeCopyNoAlloc(grid, grid_tmp);
  
  //initialize new node counter
  vtkIdType l_newNodeId = l_N_points;
  
  //actor
  for(int i=0;i<cells.size();i++) {
    if(l_marked_cells[i]==1) {
      stencil_t S = l_StencilVector[i];
      
      //calculate midpoint
      vec3_t A,B;
      grid_tmp->GetPoint(S.p[1],A.data());
      grid_tmp->GetPoint(S.p[3],B.data());
      vec3_t M=0.5*(A+B);
      
      //project point
      M=project(M);
      //add point
      grid_tmp->GetPoints()->SetPoint(l_newNodeId, M.data());
      copyNodeData(grid_tmp,S.p[1],grid_tmp,l_newNodeId);
      
      // inserted edge point = type of the edge on which it is inserted
      EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");
      node_type->SetValue(l_newNodeId, getNewNodeType(S) );
      
      if(S.twocells && S.neighbour_type==VTK_TRIANGLE){//2 triangles
        //four new triangles
        vtkIdType pts_triangle[4][3];
        for(int i=0;i<4;i++)
        {
          pts_triangle[i][0]=S.p[i];
          pts_triangle[i][1]=S.p[(i+1)%4];
          pts_triangle[i][2]=l_newNodeId;
        }
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
        grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
        
        vtkIdType newCellId;
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
        copyCellData(grid_tmp,S.id_cell2,grid_tmp,newCellId);
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
        copyCellData(grid_tmp,S.id_cell1,grid_tmp,newCellId);
      }
      else if(!S.twocells) {//1 triangle
        //two new triangles
        vtkIdType pts_triangle[2][3];
        pts_triangle[0][0]=S.p[0];
        pts_triangle[0][1]=S.p[1];
        pts_triangle[0][2]=l_newNodeId;
        pts_triangle[1][0]=S.p[3];
        pts_triangle[1][1]=S.p[0];
        pts_triangle[1][2]=l_newNodeId;
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
        
        vtkIdType newCellId;
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[1]);
        copyCellData(grid_tmp,S.id_cell1,grid_tmp,newCellId);
      }
      else {
        cout<<"I DON'T KNOW HOW TO SPLIT THIS CELL!!!"<<endl;
        EG_BUG;
      }
      
      //increment ID
      l_newNodeId++;
    }
  }
  
  //update grid
  makeCopy(grid_tmp,grid);
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_EP_all END==="<<endl;
  return(0);
}

char InsertPoints::getNewNodeType(stencil_t S)
{
  vtkIdType id_node1 = S.p[1];
  vtkIdType id_node2 = S.p[3];
  
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  if( node_type->GetValue(id_node1)==VTK_SIMPLE_VERTEX || node_type->GetValue(id_node2)==VTK_SIMPLE_VERTEX ) {
    return VTK_SIMPLE_VERTEX;
  }
  else {
    QVector <vtkIdType> PSP = getPotentialSnapPoints(id_node1);
    if( PSP.contains(id_node2) ) {
      EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
      if( cell_code->GetValue(S.id_cell1) != cell_code->GetValue(S.id_cell2) ) {
        return VTK_BOUNDARY_EDGE_VERTEX;
      }
      else {
        return VTK_FEATURE_EDGE_VERTEX;
      }
    }
    else {
      return VTK_SIMPLE_VERTEX;
    }
  }
}

///@@@ TODO:
       //============================================
      ///@@@  TODO: PRIORITY 1: Update node info (densities+type)
// EG_VTKDCN(vtkIntArray, node_specified_density, grid_tmp, "node_specified_density");//density index from table
// EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid_tmp, "node_meshdensity_desired");//what we want
// EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid_tmp, "node_meshdensity_current");//what we have
// EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");//node type
      //============================================

//       //part 1
//       node_type->SetValue(l_newNodeId,VTK_SIMPLE_VERTEX);
// 
//       //part 2
//       double total_dist=0;
//       double avg_dist=0;
//       for(int i=0;i<N_pts;i++)
//       {
//         double dist=(corner[i]-C).abs();
//         total_dist+=dist;
//         node_meshdensity_current->SetValue(pts[i],NewCurrentMeshDensity(pts[i],dist));
//       }
//       avg_dist=total_dist/(double)N_pts;
//       node_meshdensity_current->SetValue(l_newNodeId,1./avg_dist);
// 
//       //part 3
//       VertexMeshDensity nodeVMD;
//       nodeVMD.type=node_type->GetValue(l_newNodeId);
//       nodeVMD.density=0;
//       nodeVMD.CurrentNode=l_newNodeId;
//       EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
//       nodeVMD.BCmap[cell_code->GetValue(id_cell)]=2;
// 
//       int idx=VMDvector.indexOf(nodeVMD);
//       node_specified_density->SetValue(l_newNodeId, idx);
// 
//       //part 4
//       if(idx!=-1)//specified
//       {
//         node_meshdensity_desired->SetValue(l_newNodeId, VMDvector[idx].density);
//       }
//       else//unspecified
//       {
//         double D=DesiredMeshDensity(l_newNodeId);
//         node_meshdensity_desired->SetValue(l_newNodeId, D);
//       }

      //============================================
