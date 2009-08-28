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

#include "checksurfaceintegrity.h"

using namespace GeometryTools;

SwapTriangles::SwapTriangles() : SurfaceOperation()
{
  m_RespectBC = false;
  m_FeatureSwap = false;
  m_FeatureAngle = GeometryTools::deg2rad(30);
  m_MaxNumLoops = 20;
}

void SwapTriangles::operate()
{
  static int nStatic_SwapTriangles;    // Value of nStatic_SwapTriangles is retained between each function call
  nStatic_SwapTriangles++;
  
  int N_swaps;
  int N_total = 0;
  int loop = 1;
  do {
    N_swaps = 0;
    setAllSurfaceCells();
    l2g_t  cells = getPartCells();
    g2l_t _cells = getPartLocalCells();
    l2l_t  n2c   = getPartN2C();
    g2l_t _nodes = getPartLocalNodes();
    EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
    QVector<bool> l_marked(cells.size());
    foreach (vtkIdType id_cell, cells) {
      if (!m_BoundaryCodes.contains(cell_code->GetValue(id_cell)) && grid->GetCellType(id_cell) == VTK_TRIANGLE) { //if it is a selected triangle
        if (!l_marked[_cells[id_cell]]) {
          for (int j = 0; j < 3; ++j) {
            bool swap = false;
            stencil_t S = getStencil(id_cell, j);
            if(S.twocells && S.sameBC && (S.neighbour_type == VTK_TRIANGLE)) {
              if( ! isEdge(S.p[0],S.p[2]) ) {
                if (!l_marked[_cells[S.id_cell2]]) {
                  vec3_t x3[4], x3_0(0,0,0);
                  vec2_t x[4];
                  for (int k = 0; k < 4; ++k) {
                    grid->GetPoints()->GetPoint(S.p[k], x3[k].data());
                    x3_0 += x3[k];
                  }
                  vec3_t n1 = triNormal(x3[0], x3[1], x3[3]);
                  vec3_t n2 = triNormal(x3[1], x3[2], x3[3]);
                  if (m_FeatureSwap || GeometryTools::angle(n1, n2) < m_FeatureAngle) {
                    if(testSwap(S)) {
                      vec3_t n = n1 + n2;
                      n.normalise();
                      vec3_t ex = orthogonalVector(n);
                      vec3_t ey = ex.cross(n);
                      for (int k = 0; k < 4; ++k) {
                        x[k] = vec2_t(x3[k]*ex, x3[k]*ey);
                      }
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
                        }
                      } else {
                        ok = false;
                        swap = true;
                      }
                      if (ok) {
                        if ((xm1 - x[2]).abs() < (xm1 - x[0]).abs()) {
                          swap = true;
                        }
                        if ((xm2 - x[0]).abs() < (xm2 - x[2]).abs()) {
                          swap = true;
                        }
                      }
                    }// end of testswap
                  } //end of if feature angle
                } //end of if l_marked
              } //end of if TestSwap
            } //end of S valid
            
            if (swap) {
              l_marked[_cells[S.id_cell1]] = true;
              l_marked[_cells[S.id_cell2]] = true;
              for(int i = 0; i < 4; i++) {
                foreach(int i_neighbour, n2c[_nodes[S.p[i]]]) {
                  l_marked[i_neighbour] = true;
                }
              }
              vtkIdType new_pts1[3], new_pts2[3];
              new_pts1[0] = S.p[1];
              new_pts1[1] = S.p[2];
              new_pts1[2] = S.p[0];
              new_pts2[0] = S.p[2];
              new_pts2[1] = S.p[3];
              new_pts2[2] = S.p[0];
              vtkIdType old_pts1[3], old_pts2[3];
              old_pts1[0] = S.p[0];
              old_pts1[1] = S.p[1];
              old_pts1[2] = S.p[3];
              old_pts2[0] = S.p[2];
              old_pts2[1] = S.p[3];
              old_pts2[2] = S.p[1];
              grid->ReplaceCell(S.id_cell1, 3, new_pts1);
              grid->ReplaceCell(S.id_cell2, 3, new_pts2);
              ++N_swaps;
              ++N_total;
              break;
            } //end of if swap
          } //end of loop through sides
        } //end of if marked
      } //end of if selected triangle
    } //end of loop through cells
    ++loop;
  } while ((N_swaps > 0) && (loop <= m_MaxNumLoops));
  //cout << N_total << " triangles have been swapped" << endl;
}

bool SwapTriangles::testSwap(stencil_t S)
{
  // old triangles
  vec3_t n1_old = triNormal(grid, S.p[0], S.p[1], S.p[3]);
  vec3_t n2_old = triNormal(grid, S.p[2], S.p[3], S.p[1]);
  
  // new triangles
  vec3_t n1_new = triNormal(grid, S.p[1], S.p[2], S.p[0]);
  vec3_t n2_new = triNormal(grid, S.p[3], S.p[0], S.p[2]);
  
  // top point
  vec3_t x_summit(0,0,0);
  vec3_t x[4];
  double l_max = 0;
  for (int k = 0; k < 4; ++k) {
    grid->GetPoints()->GetPoint(S.p[k], x[k].data());
    x_summit += x[k];
  }
  for (int k = 0; k < 4; ++k) {
    int i = k;
    int j = k + 1;
    if (j == 4) {
      j = 0;
    }
    l_max = max(l_max, (x[i]-x[j]).abs());
  }
  x_summit *= 0.25;
  vec3_t n = n1_old + n2_old;
  n.normalise();
  x_summit += 3*l_max*n;

  // old volumes
  double V1_old = tetraVol(x[0], x[1], x[3], x_summit, true);
  double V2_old = tetraVol(x[2], x[3], x[1], x_summit, true);
  // new volumes
  double V1_new = tetraVol(x[1], x[2], x[0], x_summit, true);
  double V2_new = tetraVol(x[3], x[0], x[2], x_summit, true);

  return(V1_old>0 && V2_old>0 && V1_new>0 && V2_new>0 );
  //return(V1_new > 0 && V2_new > 0);
}

bool SwapTriangles::isEdge(vtkIdType id_node1, vtkIdType id_node2)
{
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();

  bool ret = false;
  foreach(int i_node, n2n[_nodes[id_node1]]) {
    vtkIdType id_node = nodes[i_node];
    if( id_node == id_node2 ) ret = true;
  }
  return(ret);
}
