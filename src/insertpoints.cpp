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

InsertPoints::InsertPoints() : SurfaceOperation()
{
  setQuickSave(true);
  getSet("surface meshing", "point insertion threshold", 1, m_Threshold);
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
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2c   = getPartN2C();
  l2l_t  c2c = getPartC2C();
  
  UpdatePotentialSnapPoints(true);
   
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, cl, grid, "node_meshdensity_desired");

  int num_newpoints=0;
  int num_newcells=0;
  
  QVector <int> marked_cells(cells.size(), 0);
  QVector <stencil_t> stencil_vector(cells.size());
  
  // find potential edges for splitting
  QList<edge_t> edges;
  for (int i = 0; i < cells.size(); ++i) {
    vtkIdType id_cell = cells[i];

    // if cell is selected and a triangle
    if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell)) && (grid->GetCellType(id_cell) == VTK_TRIANGLE)) {
      int j_split = -1;
      double L_max = 0;
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      //find best side to split (longest)
      for (int j = 0; j < 3; ++j) {
        //check if neighbour cell on this side is also selected
        int i_cell_neighbour = c2c[_cells[id_cell]][j];
        vtkIdType id_cell_neighbour = cells[i_cell_neighbour];
        if(m_BoundaryCodes.contains(cell_code->GetValue(id_cell_neighbour))) {
          vtkIdType id_node1 = pts[j];
          vtkIdType id_node2 = pts[(j+1)%N_pts];
          double L  = distance(grid, id_node1, id_node2);
          double L1 = cl->GetValue(id_node1);
          double L2 = cl->GetValue(id_node2);
          if (L > m_Threshold*min(L1,L2)) {
            if (L > L_max) {
              j_split = j;
              L_max = L;
            }
          }
        }
      }
      if (j_split != -1) {
        stencil_t S = getStencil(id_cell, j_split);
        edge_t E;
        E.S = S;
        E.L1 = cl->GetValue(S.p[1]);
        E.L2 = cl->GetValue(S.p[3]);
        E.L12 = distance(grid, S.p[1], S.p[3]);
        edges.push_back(E);
      }
    }
  }

  qSort(edges);

  //counter
  foreach (edge_t E, edges) {
    if (E.S.twocells && (E.S.neighbour_type == VTK_TRIANGLE)) {
      int i_cells1 = _cells[E.S.id_cell1];
      int i_cells2 = _cells[E.S.id_cell2];
      if (!marked_cells[i_cells1]  && !marked_cells[i_cells2]) {
        stencil_vector[i_cells1] = E.S;
        marked_cells[i_cells1] = 1;
        marked_cells[i_cells2] = 2;
        ++num_newpoints;
        num_newcells += 2;
      }
    } else if (!E.S.twocells) {
      int i_cells1 = _cells[E.S.id_cell1];
      if (!marked_cells[i_cells1]) {
        stencil_vector[i_cells1] = E.S;
        marked_cells[i_cells1] = 1;
        ++num_newpoints;
        num_newcells += 1;
      }
    }
  }

  //initialize grid_tmp
  int l_N_points=grid->GetNumberOfPoints();
  int l_N_cells=grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp, l_N_cells + num_newcells, l_N_points + num_newpoints);
  makeCopyNoAlloc(grid, grid_tmp);
  
  //initialize new node counter
  vtkIdType id_new_node = l_N_points;
  
  //actor
  for (int i = 0; i < cells.size(); i++) {
    if (marked_cells[i] == 1) {
      stencil_t S = stencil_vector[i];
      
      //calculate midpoint
      vec3_t A,B;
      grid_tmp->GetPoint(S.p[1],A.data());
      grid_tmp->GetPoint(S.p[3],B.data());
      vec3_t M=0.5*(A+B);
      
      //add point
      grid_tmp->GetPoints()->SetPoint(id_new_node, M.data());
      copyNodeData(grid_tmp,S.p[1],grid_tmp,id_new_node); ///@@@ maybe trouble
      
      // inserted edge point = type of the edge on which it is inserted
      EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");
      node_type->SetValue(id_new_node, getNewNodeType(S) );
      
      if(S.twocells && S.neighbour_type==VTK_TRIANGLE) { //2 triangles
        //four new triangles
        vtkIdType pts_triangle[4][3];
        for(int i = 0; i < 4; ++i) {
          pts_triangle[i][0] = S.p[i];
          pts_triangle[i][1] = S.p[(i+1)%4];
          pts_triangle[i][2] = id_new_node;
        }
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
        grid_tmp->ReplaceCell(S.id_cell2 , 3, pts_triangle[1]);
        
        vtkIdType newCellId;
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[2]);
        copyCellData(grid_tmp,S.id_cell2,grid_tmp,newCellId);
        
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[3]);
        copyCellData(grid_tmp,S.id_cell1,grid_tmp,newCellId);
      } else if(!S.twocells) { //1 triangle
        //two new triangles
        vtkIdType pts_triangle[2][3];
        pts_triangle[0][0] = S.p[0];
        pts_triangle[0][1] = S.p[1];
        pts_triangle[0][2] = id_new_node;
        pts_triangle[1][0] = S.p[3];
        pts_triangle[1][1] = S.p[0];
        pts_triangle[1][2] = id_new_node;
        
        grid_tmp->ReplaceCell(S.id_cell1 , 3, pts_triangle[0]);
        
        vtkIdType newCellId;
        newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[1]);
        copyCellData(grid_tmp,S.id_cell1,grid_tmp,newCellId);
      } else {
        cout<<"I DON'T KNOW HOW TO SPLIT THIS CELL!!!"<<endl;
        EG_BUG;
      }
      
      //increment ID
      id_new_node++;
    }
  }
  
  //update grid
  makeCopy(grid_tmp,grid);
  
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

