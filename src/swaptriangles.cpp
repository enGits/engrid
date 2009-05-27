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
#include "swaptriangles.h"
#include "guimainwindow.h"

#include <vtkCellArray.h>

using namespace GeometryTools;

SwapTriangles::SwapTriangles()
: Operation()
{
  m_RespectBC=false;
  m_FeatureSwap=false;
}

void SwapTriangles::prepare()
{
  getAllCellsOfType(VTK_TRIANGLE, cells, grid);
  QList<vtkIdType> ex_cells;
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  foreach (vtkIdType id_cell, cells) {
    if (!boundary_codes.contains(bc->GetValue(id_cell))) {
      ex_cells.append(id_cell);
    };
    if (grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      EG_BUG;
    };
  };
  cells.resize(ex_cells.size());
  qCopy(ex_cells.begin(), ex_cells.end(), cells.begin());
  createCellMapping(cells, _cells, grid);
  createCellToCell(cells, c2c, grid);
  marked.resize(cells.size());
};

void SwapTriangles::operate()
{
  cout << "swapping edges of boundary triangles (Delaunay)" << endl;
  using namespace GeometryTools;
  prepare();
  int N_swaps;
  int N_total = 0;
  int loop = 1;
  do {
    N_swaps = 0;
    createCellToCell(cells, c2c, grid);
    //NOTE: This for loop can eventually be removed because if undefined, it's probably false.
    for (int i_cells = 0; i_cells < cells.size(); ++i_cells) {
      marked[i_cells] = false;
    };
    foreach (vtkIdType id_cell, cells) {
//       cout<<"===>id_cell="<<id_cell<<" loop="<<loop<<endl;
      if (!marked[_cells[id_cell]]) {
        for (int j = 0; j < 3; ++j) {
          bool swap = false;
          stencil_t S = getStencil(id_cell, j, m_RespectBC);
          if (S.valid) {
            if (!marked[_cells[S.id_cell2]]) {
              vec3_t x3[4], x3_0(0,0,0);
              vec2_t x[4];
              for (int k = 0; k < 4; ++k) {
                grid->GetPoints()->GetPoint(S.p[k], x3[k].data());
                x3_0 += x3[k];
              };
              vec3_t n1 = triNormal(x3[0], x3[1], x3[3]);
              vec3_t n2 = triNormal(x3[1], x3[2], x3[3]);
/*              n1.normalise();
              n2.normalise();*/
              if ( m_FeatureSwap || (n1*n2) > 0.8*n1.abs()*n2.abs() ) {
                vec3_t n = n1 + n2;
                n.normalise();
                vec3_t ex = orthogonalVector(n);
                vec3_t ey = ex.cross(n);
                for (int k = 0; k < 4; ++k) {
                  x[k] = vec2_t(x3[k]*ex, x3[k]*ey);
                };
                vec2_t r1, r2, r3, u1, u2, u3;
                r1 = 0.5*(x[0] + x[1]); u1 = turnLeft(x[1] - x[0]);
                r2 = 0.5*(x[1] + x[2]); u2 = turnLeft(x[2] - x[1]);
                r3 = 0.5*(x[1] + x[3]); u3 = turnLeft(x[3] - x[1]);
                double k, l;
                vec2_t xm1, xm2;
                bool ok = true;
                if (intersection(k, l, r1, u1, r3, u3)) {
                  xm1 = r1 + k*u1;
                  if (intersection(k, l, r2, u2, r3, u3)) {
                    xm2 = r2 + k*u2;
                  } else {
                    ok = false;
                  };
                } else {
                  ok = false;
                  swap = true;
                };
                if (ok) {
                  if ((xm1 - x[2]).abs() < (xm1 - x[0]).abs()) {
                    swap = true;
                  };
                  if ((xm2 - x[0]).abs() < (xm2 - x[2]).abs()) {
                    swap = true;
                  };
                };
              };
            };
          };
          if (swap) {
            marked[_cells[S.id_cell1]] = true;
            marked[_cells[S.id_cell2]] = true;
            if (c2c[_cells[S.id_cell1]][0] != -1) {
              marked[c2c[_cells[S.id_cell1]][0]] = true;
            } else if (c2c[_cells[S.id_cell1]][1] != -1) {
              marked[c2c[_cells[S.id_cell1]][1]] = true;
            } else if (c2c[_cells[S.id_cell1]][2] != -1) {
              marked[c2c[_cells[S.id_cell1]][2]] = true;
            };
            if (c2c[_cells[S.id_cell2]][0] != -1) {
              marked[c2c[_cells[S.id_cell2]][0]] = true;
            } else if (c2c[_cells[S.id_cell2]][1] != -1) {
              marked[c2c[_cells[S.id_cell2]][1]] = true;
            } else if (c2c[_cells[S.id_cell2]][2] != -1) {
              marked[c2c[_cells[S.id_cell2]][2]] = true;
            };
            vtkIdType new_pts1[3], new_pts2[3];
            new_pts1[0] = S.p[1];
            new_pts1[1] = S.p[2];
            new_pts1[2] = S.p[0];
            new_pts2[0] = S.p[2];
            new_pts2[1] = S.p[3];
            new_pts2[2] = S.p[0];
            grid->ReplaceCell(S.id_cell1, 3, new_pts1);
            grid->ReplaceCell(S.id_cell2, 3, new_pts2);
            ++N_swaps;
            ++N_total;
//             GuiMainWindow::pointer()->QuickSave();
/*            if(S.id_cell1==138) {
              cout<<"swapping!!!!!!"<<endl;
              cout<<"S.id_cell1="<<S.id_cell1<<endl;
              cout<<"S.id_cell2="<<S.id_cell2<<endl;
              cout<<"marked[_cells[S.id_cell1]]="<<marked[_cells[S.id_cell1]]<<endl;
              cout<<"marked[_cells[S.id_cell2]]="<<marked[_cells[S.id_cell2]]<<endl;
            }*/
            break;
          };
        };
      };
    };
    ++loop;
  } while ((N_swaps > 0) && (loop <= 20));
  cout << N_total << " triangles have been swapped" << endl;
};

bool SwapTriangles::TestSwap(stencil_t S)
{
  //old triangles
  vec3_t n1_old=triNormal(grid,S.p[0],S.p[1],S.p[3]);
  vec3_t n2_old=triNormal(grid,S.p[2],S.p[3],S.p[1]);
  
  //new triangles
  vec3_t n1_new=triNormal(grid,S.p[1],S.p[2],S.p[0]);
  vec3_t n2_new=triNormal(grid,S.p[3],S.p[0],S.p[2]);
  
  //top point
  vec3_t Summit=n1_old+n2_old;
  vec3_t M[4];
  for (int k = 0; k < 4; ++k) {
    grid->GetPoints()->GetPoint(S.p[k], M[k].data());
  };
  
  //old volumes
  double V1_old=tetraVol(M[0], Summit, M[1], M[3], true);
  double V2_old=tetraVol(M[2], Summit, M[3], M[1], true);
  //new volumes
  double V1_new=tetraVol(M[1], Summit, M[2], M[0], true);
  double V2_new=tetraVol(M[3], Summit, M[0], M[2], true);

  
/*  if(S.id_cell1==138)
  {
    cout<<"V1_old="<<V1_old<<endl;
    cout<<"V2_old="<<V2_old<<endl;
    cout<<"V1_new="<<V1_new<<endl;
    cout<<"V2_new="<<V2_new<<endl;
    
    cout<<"n1_old="<<n1_old<<endl;
    cout<<"n2_old="<<n2_old<<endl;
    cout<<"n1_new="<<n1_new<<endl;
    cout<<"n2_new="<<n2_new<<endl;
    
    cout<<"n1_old*n1_old="<<n1_old*n1_old<<endl;
    cout<<"n2_old*n1_old="<<n2_old*n1_old<<endl;
    cout<<"n1_old*n2_old="<<n1_old*n2_old<<endl;
    cout<<"n2_old*n2_old="<<n2_old*n2_old<<endl;
    
    cout<<"n1_old*n1_new="<<n1_old*n1_new<<endl;
    cout<<"n2_old*n1_new="<<n2_old*n1_new<<endl;
    cout<<"n1_old*n2_new="<<n1_old*n2_new<<endl;
    cout<<"n2_old*n2_new="<<n2_old*n2_new<<endl;
    
    cout<<"n1_new*n1_old="<<n1_new*n1_old<<endl;
    cout<<"n2_new*n1_old="<<n2_new*n1_old<<endl;
    cout<<"n1_new*n2_old="<<n1_new*n2_old<<endl;
    cout<<"n2_new*n2_old="<<n2_new*n2_old<<endl;
    
    cout<<"n1_new*n1_new="<<n1_new*n1_new<<endl;
    cout<<"n2_new*n1_new="<<n2_new*n1_new<<endl;
    cout<<"n1_new*n2_new="<<n1_new*n2_new<<endl;
    cout<<"n2_new*n2_new="<<n2_new*n2_new<<endl;
  }*/
//   return(n1_new*n1_old>0 && n1_new*n2_old>0 && n2_new*n1_old>0 && n2_new*n2_old>0);
  return(V1_old>0 && V2_old>0 && V1_new>0 && V2_new>0 );
}
