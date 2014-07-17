// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2014 enGits GmbH                                      +
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
  int N1 = m_Grid->GetNumberOfPoints();
  insertPoints();
  int N2 = m_Grid->GetNumberOfPoints();
  m_NumInserted = N2 - N1;
}

///\todo Adapt this code for multiple volumes.
int InsertPoints::insertPoints()
{
  setAllSurfaceCells();
  l2g_t  cells = getPartCells();
  g2l_t _cells = getPartLocalCells();

  updateNodeInfo();

  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  EG_VTKDCN(vtkDoubleArray, characteristic_length_desired, m_Grid, "node_meshdensity_desired");

  int num_newpoints=0;
  int num_newcells=0;

  QVector <int> marked_cells(cells.size(), 0);
  QVector <stencil_t> stencil_vector(cells.size());

  // find potential edges for splitting
  QList<edge_t> edges;
  for (int i = 0; i < cells.size(); ++i) {
    vtkIdType id_cell = cells[i];

    // if cell is selected and a triangle
    if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell)) && (m_Grid->GetCellType(id_cell) == VTK_TRIANGLE)) {
      int j_split = -1;
      double L_max = 0;
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);

      //find best side to split (longest)
      for (int j = 0; j < 3; ++j) {

        //check if neighbour cell on this side is also selected
        stencil_t S = getStencil(id_cell, j);
        bool selected_edge = true;
        for (int i_cell_neighbour=1;i_cell_neighbour<S.id_cell.size();i_cell_neighbour++) {
          vtkIdType id_cell_neighbour = S.id_cell[i_cell_neighbour];
          if (!m_BoundaryCodes.contains(cell_code->GetValue(id_cell_neighbour)) || S.type_cell[i_cell_neighbour] != VTK_TRIANGLE) {
            selected_edge=false;
          }
        } // end of loop through neighbour cells

        if (selected_edge) {
          vtkIdType id_node1 = pts[j];
          vtkIdType id_node2 = pts[(j+1)%N_pts];
          double L  = distance(m_Grid, id_node1, id_node2);
          double L1 = characteristic_length_desired->GetValue(id_node1);
          double L2 = characteristic_length_desired->GetValue(id_node2);
          if (L > m_Threshold*min(L1,L2)) {
            if (L > L_max) {
              j_split = j;
              L_max = L;
            }
          }
        }
      } // end of loop through sides
      if (j_split != -1) {
        stencil_t S = getStencil(id_cell, j_split);
        edge_t E;
        E.S = S;
        E.L1 = characteristic_length_desired->GetValue(S.p1);
        E.L2 = characteristic_length_desired->GetValue(S.p2);
        E.L12 = distance(m_Grid, S.p1, S.p2);
        edges.push_back(E);
      }
    }
  }

  qSort(edges);

  //counter
  foreach (edge_t E, edges) {
    int i_cells1 = _cells[E.S.id_cell[0]];
    bool all_unmarked = true;
    for (int i = 0; i < E.S.id_cell.size(); i++) {
      int i_cells = _cells[E.S.id_cell[i]];
      if (marked_cells[i_cells]) {
        all_unmarked=false;
      }
    }
    if (all_unmarked) {
      stencil_vector[i_cells1] = E.S;
      marked_cells[i_cells1] = 1;
      for (int i = 1; i < E.S.id_cell.size(); i++) {
        int i_cells = _cells[E.S.id_cell[i]];
        marked_cells[i_cells] = 2;
      }
      ++num_newpoints;
      num_newcells += E.S.id_cell.size();
    }
  }

  //initialize grid_tmp
  int l_N_points=m_Grid->GetNumberOfPoints();
  int l_N_cells=m_Grid->GetNumberOfCells();
  EG_VTKSP(vtkUnstructuredGrid,grid_tmp);
  allocateGrid(grid_tmp, l_N_cells + num_newcells, l_N_points + num_newpoints);
  makeCopyNoAlloc(m_Grid, grid_tmp);

  //initialize new node counter
  vtkIdType id_new_node = l_N_points;

  //actor
  for (int i = 0; i < cells.size(); i++) {
    if (marked_cells[i] == 1) {
      stencil_t S = stencil_vector[i];

      //calculate midpoint
      vec3_t A,B;
      grid_tmp->GetPoint(S.p1,A.data());
      grid_tmp->GetPoint(S.p2,B.data());
      vec3_t M=0.5*(A+B);

      //add point
      grid_tmp->GetPoints()->SetPoint(id_new_node, M.data());
      copyNodeData(grid_tmp,S.p1,grid_tmp,id_new_node); ///\todo maybe trouble

      // inserted edge point = type of the edge on which it is inserted
      EG_VTKDCN(vtkCharArray, node_type, grid_tmp, "node_type");
      node_type->SetValue(id_new_node, getNewNodeType(S) );

      // insert new cells
      //four new triangles
      int N = S.id_cell.size();
      QVector< QVector<vtkIdType> > pts_triangle(2*N,QVector<vtkIdType>(3));

      for(int i_triangle = 0; i_triangle < N; i_triangle++) {
        vtkIdType *pts, N_pts;
        grid_tmp->GetCellPoints(S.id_cell[i_triangle], N_pts, pts);

        bool direct;
        for(int i_pts = 0; i_pts < N_pts; i_pts++) {
          if (pts[i_pts] == S.p1 ) {
            if (pts[(i_pts+1)%N_pts] == S.p2 ) {
              direct = true;
            } else {
              direct = false;
            }
          }
        }
        if (direct) {
          pts_triangle[i_triangle][0] = S.p1;
          pts_triangle[i_triangle][1] = id_new_node;
          pts_triangle[i_triangle][2] = S.id_node[i_triangle];

          pts_triangle[i_triangle+N][0] = id_new_node;
          pts_triangle[i_triangle+N][1] = S.p2;
          pts_triangle[i_triangle+N][2] = S.id_node[i_triangle];
        } else {
          pts_triangle[i_triangle][0] = S.p2;
          pts_triangle[i_triangle][1] = id_new_node;
          pts_triangle[i_triangle][2] = S.id_node[i_triangle];

          pts_triangle[i_triangle+N][0] = id_new_node;
          pts_triangle[i_triangle+N][1] = S.p1;
          pts_triangle[i_triangle+N][2] = S.id_node[i_triangle];
        }

        grid_tmp->ReplaceCell(S.id_cell[i_triangle] , 3, pts_triangle[i_triangle].data());
        vtkIdType newCellId = grid_tmp->InsertNextCell(VTK_TRIANGLE,3,pts_triangle[i_triangle+N].data());
        copyCellData(grid_tmp,S.id_cell[i_triangle],grid_tmp,newCellId);
      }

      //increment ID
      id_new_node++;
    }
  }

  //update grid
  makeCopy(grid_tmp,m_Grid);

  return(0);
}

char InsertPoints::getNewNodeType(stencil_t S)
{
  vtkIdType id_node1 = S.p1;
  vtkIdType id_node2 = S.p2;

  EG_VTKDCN(vtkCharArray, node_type, m_Grid, "node_type");
  char type1 = node_type->GetValue(id_node1);
  char type2 = node_type->GetValue(id_node2);
  if (type1 == EG_SIMPLE_VERTEX || type2 == EG_SIMPLE_VERTEX) {
    return EG_SIMPLE_VERTEX;
  } else {
    QVector <vtkIdType> PSP = getPotentialSnapPoints(id_node1);
    if (PSP.contains(id_node2)) {
      EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
      if (S.id_cell.size() < 1) {
        return EG_BOUNDARY_EDGE_VERTEX;
      } else if (S.id_cell.size() == 1) {
        EG_ERR_RETURN("Invalid surface mesh. Check this with 'Tools -> Check surface integrity'.");
        return EG_FEATURE_EDGE_VERTEX; //at best, this would be a feature edge, since it's loose.
      } else {
        if (cell_code->GetValue(S.id_cell[0]) != cell_code->GetValue(S.id_cell[1])) {
          return EG_BOUNDARY_EDGE_VERTEX;
        } else {
          if (isFeatureNode(id_node1) || isFeatureNode(id_node2)) {
            return EG_FEATURE_EDGE_VERTEX;
          } else {
            return EG_SIMPLE_VERTEX;
          }
        }
      }
    } else {
      return EG_SIMPLE_VERTEX;
    }
  }
}

