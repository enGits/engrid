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
  int N1 = grid->GetNumberOfPoints();
  insertPoints();
  int N2 = grid->GetNumberOfPoints();
  m_NumInserted = N2 - N1;
}

int InsertPoints::insertPoints()
{
  QTime start = QTime::currentTime();
  
  setAllSurfaceCells();
  l2g_t  cells = getPartCells();
  g2l_t _cells = getPartLocalCells();
  
  UpdatePotentialSnapPoints(true);
  
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  
  int num_newpoints=0;
  int num_newcells=0;
  
  QVector <int> marked_cells(cells.size(), 0);
  QVector <stencil_t> stencil_vector(cells.size());
  
  //counter
  for (int i = 0; i < cells.size(); ++i) {
    vtkIdType id_cell = cells[i];
    if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell)) && (grid->GetCellType(id_cell) == VTK_TRIANGLE)) {//if selected and triangle cell
      int j_split = -1;
      double L_max = 0;
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int j = 0; j < 3; ++j) {
        vtkIdType id_node1 = pts[j];
        vtkIdType id_node2 = pts[(j+1)%N_pts];
        double L = distance(grid, id_node1, id_node2);
        if (L > max(desiredEdgeLength(id_node1), desiredEdgeLength(id_node2))) {
          if (L > L_max) {
            j_split = j;
            L_max = L;
          }
        }
      }
      if (j_split != -1) {
        stencil_t S = getStencil(id_cell, j_split);
        if (S.twocells && (S.neighbour_type == VTK_TRIANGLE)) {
          if (!marked_cells[i]  && !marked_cells[_cells[S.id_cell2]]) {
            stencil_vector[i] = S;
            marked_cells[i] = 1;
            marked_cells[_cells[S.id_cell2]] = 2;
            num_newpoints++;
            num_newcells += 2;
          }
        } else if (!S.twocells) {
          if (!marked_cells[i]) {
            stencil_vector[i] = S;
            marked_cells[i] = 1;
            num_newpoints++;
            num_newcells+=1;
          }
        }
      } //end of loop through sides
    } //end of if selected and triangle cell
  } //end of counter loop
  
  //initialize grid_tmp
  int l_N_points=grid->GetNumberOfPoints();
  int l_N_cells=grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp, l_N_cells + num_newcells, l_N_points + num_newpoints);
  makeCopyNoAlloc(grid, grid_tmp);
  
  //initialize new node counter
  vtkIdType l_newNodeId = l_N_points;
  
  //actor
  for (int i = 0; i < cells.size(); i++) {
    if (marked_cells[i] == 1) {
      stencil_t S = stencil_vector[i];
      
      //calculate midpoint
      vec3_t A,B;
      grid_tmp->GetPoint(S.p[1],A.data());
      grid_tmp->GetPoint(S.p[3],B.data());
      vec3_t M=0.5*(A+B);
      
      //project point
      //M=project(M);
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
  
  //cout << start.msecsTo(QTime::currentTime()) << " milliseconds elapsed" << endl;
  //cout<<"===insert_EP_all END==="<<endl;
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
      } else {
        return VTK_FEATURE_EDGE_VERTEX;
      }
    } else {
      return VTK_SIMPLE_VERTEX;
    }
  }
}

