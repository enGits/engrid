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

#include <vtkCharArray.h>

#include <QTime>

InsertPoints::InsertPoints()
 : Operation()
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
  bool result=L_k(j,K)>0.5*(G_k(j)+G_k(K));
  if(DebugLevel>0 && result){
    cout<<"j="<<j<<endl;
    cout<<"K="<<K<<endl;
    cout<<"G_k(j)="<<G_k(j)<<endl;
    cout<<"G_k(K)="<<G_k(K)<<endl;
    cout<<"0.5*(G_k(j)+G_k(K))="<<0.5*(G_k(j)+G_k(K))<<endl;
    cout<<"L_k(j,K)="<<L_k(j,K)<<endl;
  }
  return ( result );
}

bool InsertPoints::SplitSide(vtkIdType id_cell,int side)
{
  vtkIdType N_pts,*pts;
  grid->GetCellPoints(id_cell,N_pts,pts);
  return( insert_edgepoint(pts[side],pts[(side+1)%N_pts]) );
}

int InsertPoints::insert_FP_counter()
{
  cout<<"===insert_FP_counter() START==="<<endl;
  QTime start = QTime::currentTime();
  
  int l_N_inserted_FP=0;

  ///@@@  TODO: optimize
  //unmark cells and nodes
  m_marked_cells.clear();
  
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
    if( !m_marked_cells[id_cell] && insert_fieldpoint(id_cell) )
    {
      l_N_inserted_FP++;
      m_marked_cells[id_cell]=true;
      m_N_newcells+=2;
      m_N_newpoints+=1;
    }
  }
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_FP_counter() END==="<<endl;
  return(l_N_inserted_FP);
}

int InsertPoints::insert_FP_actor(vtkUnstructuredGrid* grid_tmp)
{
  cout<<"===insert_FP_actor START==="<<endl;
  QTime start = QTime::currentTime();
  
  //initialize new node counter
  vtkIdType l_newNodeId=m_N_points;
  
  vtkCellLocator* l_CellLocator = vtkCellLocator::New();
  l_CellLocator->SetDataSet(m_ProjectionSurface);
  l_CellLocator->BuildLocator();
  
  ///@@@  TODO: optimize
  //unmark cells
  m_marked_cells.clear();//why?
  
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  foreach(vtkIdType id_cell, m_SelectedCells)
  {
    if( !m_marked_cells[id_cell] && insert_fieldpoint(id_cell) )
    {
      m_marked_cells[id_cell]=true;
      
      vtkIdType newBC=cell_code_tmp->GetValue(id_cell);
      
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      vec3_t C(0,0,0);
      
      int N_neighbours=N_pts;
      if(DebugLevel>42) cout<<"N_neighbours="<<N_neighbours<<endl;
      vec3_t corner[4];
      vtkIdType pts_triangle[4][3];
      for(int i=0;i<N_neighbours;i++)
      {
        grid->GetPoints()->GetPoint(pts[i], corner[i].data());
        C+=corner[i];
      }
      C=(1/(double)N_neighbours)*C;
      
      //============================================
      // ADD POINT
      C=project(C);
      addPoint(grid_tmp,l_newNodeId,C.data(),m_CellLocator);
      

      //============================================
      ///@@@  TODO: PRIORITY 1: Update node info (densities+type)
      EG_VTKDCN(vtkIntArray, node_specified_density, grid_tmp, "node_specified_density");//density index from table
      EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid_tmp, "node_meshdensity_desired");//what we want
      EG_VTKDCN(vtkDoubleArray, node_meshdensity_current, grid_tmp, "node_meshdensity_current");//what we have
      EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");//node type
      //============================================
      
//       //part 1
//       node_type->SetValue(l_newNodeId,VTK_SIMPLE_VERTEX);
// 
//       //part 2
//       double total_dist=0;
//       double avg_dist=0;
//       for(int i=0;i<N_neighbours;i++)
//       {
//         double dist=(corner[i]-C).abs();
//         total_dist+=dist;
//         node_meshdensity_current->SetValue(pts[i],NewCurrentMeshDensity(pts[i],dist));
//       }
//       avg_dist=total_dist/(double)N_neighbours;
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

      vtkIdType intmidpoint=l_newNodeId;
      l_newNodeId++;
      
      for(int i=0;i<N_neighbours;i++)
      {
        pts_triangle[i][0]=pts[i];
        pts_triangle[i][1]=pts[(i+1)%N_neighbours];
        pts_triangle[i][2]=intmidpoint;
        if(i==0)
        {
          grid_tmp->ReplaceCell(id_cell , 3, pts_triangle[0]);
          if(cellVA(grid_tmp,id_cell)<10e-6) EG_BUG;
          cell_code_tmp->SetValue(id_cell, newBC);
        }
        else
        {
          vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[i]);
          if(cellVA(grid_tmp,newCellId)<10e-6) EG_BUG;
          cell_code_tmp->SetValue(newCellId, newBC);
        }
      }
      
    }
  }
  
  l_CellLocator->Delete();
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_FP_actor END==="<<endl;
  return(0);
}

int InsertPoints::insert_FP_all()
{
  cout<<"===insert_FP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  getAllSurfaceCells(m_AllCells,grid);
  getSurfaceCells(m_bcs, m_SelectedCells, grid);
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  getSurfaceNodes(m_bcs,m_SelectedNodes,grid);
  getNodesFromCells(m_AllCells, nodes, grid);
  setGrid(grid);
  setCells(m_AllCells);
  cout<<"m_AllCells.size()="<<m_AllCells.size()<<endl;
  
  m_N_points=grid->GetNumberOfPoints();
  m_N_cells=grid->GetNumberOfCells();
  m_N_newpoints=0;
  m_N_newcells=0;
  
  insert_FP_counter();
  
    //init grid_tmp
  m_N_points=grid->GetNumberOfPoints();
  m_N_cells=grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,m_N_cells+m_N_newcells,m_N_points+m_N_newpoints);
  m_total_N_newpoints+=m_N_newpoints; m_total_N_newcells+=m_N_newcells;
  
  makeCopyNoAlloc(grid, grid_tmp);
  
  insert_FP_actor(grid_tmp);
  
  makeCopy(grid_tmp,grid);
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_FP_all END==="<<endl;
  return(0);
}

int InsertPoints::insert_EP_all()
{
  cout<<"===insert_EP_all START==="<<endl;
  QTime start = QTime::currentTime();
  
  setAllSurfaceCells();
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, node_meshdensity_desired, grid, "node_meshdensity_desired");
  
  int l_N_newpoints=0;
  int l_N_newcells=0;
  
  QVector <int> l_marked_cells(cells.size());
  QVector <stencil_t> l_StencilVector(cells.size());
  
  //counter
  for(int i=0;i<cells.size();i++) {
    vtkIdType id_cell=cells[i];
    if(m_bcs.contains(cell_code->GetValue(id_cell))) {
      for(int j = 0; j < 3; ++j) {
        stencil_t S = getStencil(id_cell,j,false);
        if( l_marked_cells[i]==0 && l_marked_cells[_cells[S.id_cell2]]==0 ) {
          if( SplitSide(id_cell,j)) {
            l_StencilVector[i]=S;
            l_marked_cells[i]=1;
            if(cells.indexOf(S.id_cell2)>=0) {
              l_marked_cells[_cells[S.id_cell2]]=2;
            }
            l_N_newpoints++;
            if(S.valid) {
              l_N_newcells+=2;
            }
            else {
              l_N_newcells+=1;
            }
          }//end of if SplitSide
        }//end of if unmarked
      }//end of loop through sides
    }//end of if selected
  }//end of counter loop
  
  //initialize grid_tmp
  int l_N_points=grid->GetNumberOfPoints();
  int l_N_cells=grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp,l_N_cells+l_N_newcells,l_N_points+l_N_newpoints);
  cout<<"=== ALLOCATING: grid ==="<<endl;cout_grid(cout,grid);
  cout<<"=== ALLOCATING: grid_tmp ==="<<endl;cout_grid(cout,grid_tmp);
  makeCopyNoAlloc(grid, grid_tmp);
  EG_VTKDCC(vtkIntArray, cell_code_tmp, grid_tmp, "cell_code");
  
  //initialize new node counter
  vtkIdType l_newNodeId = l_N_points;
  
  //initialize vtkCellLocator
  vtkCellLocator* l_CellLocator = vtkCellLocator::New();
  l_CellLocator->SetDataSet(m_ProjectionSurface);
  l_CellLocator->BuildLocator();
  
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
      addPoint(grid_tmp,l_newNodeId,M.data(),l_CellLocator);
      
      if(S.valid){//there is a neighbour cell
        //four new triangles
        vtkIdType pts_triangle[4][3];
        for(int i=0;i<4;i++)
        {
          pts_triangle[i][0]=S.p[i];
          pts_triangle[i][1]=S.p[(i+1)%4];
          pts_triangle[i][2]=l_newNodeId;
        }
        
        int bc1=cell_code_tmp->GetValue(S.id_cell1);
        int bc2=cell_code_tmp->GetValue(S.id_cell2);
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
//         if(cellVA(grid_tmp,S.id_cell1)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(S.id_cell1, bc1);
        
        grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
//         if(cellVA(grid_tmp,S.id_cell2)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(S.id_cell2, bc2);
        
        vtkIdType newCellId;
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
//         if(cellVA(grid_tmp,newCellId)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(newCellId, bc2);
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
//         if(cellVA(grid_tmp,newCellId)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(newCellId, bc1);
      }
      else{//there is no neighbour cell
        //two new triangles
        vtkIdType pts_triangle[2][3];
        pts_triangle[0][0]=S.p[0];
        pts_triangle[0][1]=S.p[1];
        pts_triangle[0][2]=l_newNodeId;
        pts_triangle[1][0]=S.p[3];
        pts_triangle[1][1]=S.p[0];
        pts_triangle[1][2]=l_newNodeId;
        
        int bc1=cell_code_tmp->GetValue(S.id_cell1);
        cout<<"bc1="<<bc1<<endl;
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
        if(cellVA(grid_tmp,S.id_cell1)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(S.id_cell1, bc1);
        
        vtkIdType newCellId;
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[1]);
        if(cellVA(grid_tmp,newCellId)<10e-6) EG_BUG;
        cell_code_tmp->SetValue(newCellId, bc1);
      }
      
      //increment ID
      l_newNodeId++;
    }
  }
  
  //delete vtkCellLocator
  l_CellLocator->Delete();
  
  //update grid
  makeCopy(grid_tmp,grid);
  
  cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  cout<<"===insert_EP_all END==="<<endl;
  return(0);
}

double InsertPoints::NewCurrentMeshDensity(vtkIdType a_vertex,double a_dist)
{
  double total_dist=0;
  double avg_dist=0;
  int N=n2n_func(a_vertex).size();
  vec3_t C;
  grid->GetPoint(a_vertex, C.data());
  foreach(int i,n2n_func(a_vertex))
  {
    vec3_t M;
    grid->GetPoint(i, M.data());
    total_dist+=(M-C).abs();
  }
  avg_dist=(total_dist+a_dist)/(double)(N+1);
  return(1./avg_dist);
}
