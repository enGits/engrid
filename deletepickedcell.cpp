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

#include "deletepickedcell.h"
#include "guimainwindow.h"
#include "geometrytools.h"
using namespace GeometryTools;

#include <vtkCell.h>
#include <vtkGenericCell.h>

bool DeletePickedCell::SwapCell(vtkUnstructuredGrid* a_grid, stencil_t S)
{
  bool swap = false;
  if (S.valid) {
    vec3_t x3[4], x3_0(0,0,0);
    vec2_t x[4];
    for (int k = 0; k < 4; ++k) {
      a_grid->GetPoints()->GetPoint(S.p[k], x3[k].data());
      x3_0 += x3[k];
    };
    vec3_t n1 = triNormal(x3[0], x3[1], x3[3]);
    vec3_t n2 = triNormal(x3[1], x3[2], x3[3]);
    n1.normalise();
    n2.normalise();
    if ( (n1*n2) > 0.8) {
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
  if (swap) {
    vtkIdType new_pts1[3], new_pts2[3];
    new_pts1[0] = S.p[1];
    new_pts1[1] = S.p[2];
    new_pts1[2] = S.p[0];
    new_pts2[0] = S.p[2];
    new_pts2[1] = S.p[3];
    new_pts2[2] = S.p[0];
    a_grid->ReplaceCell(S.id_cell1, 3, new_pts1);
    a_grid->ReplaceCell(S.id_cell2, 3, new_pts2);
  };
  return(swap);
}

void DeletePickedCell::foobar(vtkUnstructuredGrid* src,vtkIdType quadcell)
{
  vtkIdType type_cell = grid->GetCellType(quadcell);
  cout<<"It's a "<<type_cell<<endl;
  if(type_cell==VTK_QUAD)
  {
    cout_grid(cout,src,true,true,true,true);
    EG_VTKSP(vtkUnstructuredGrid, dst);
    //src grid info
    int N_points=src->GetNumberOfPoints();
    int N_cells=src->GetNumberOfCells();
    allocateGrid(dst,N_cells+1,N_points);
    
    for (vtkIdType id_node = 0; id_node < src->GetNumberOfPoints(); ++id_node) {
      vec3_t x;
      src->GetPoints()->GetPoint(id_node, x.data());
      dst->GetPoints()->SetPoint(id_node, x.data());
      copyNodeData(src, id_node, dst, id_node);
    };
    for (vtkIdType id_cell = 0; id_cell < src->GetNumberOfCells(); ++id_cell) {
      vtkIdType N_pts, *pts;
      vtkIdType type_cell = src->GetCellType(id_cell);
      src->GetCellPoints(id_cell, N_pts, pts);
      vtkIdType id_new_cell;
      vtkIdType id_new_cell1;
      vtkIdType id_new_cell2;
      if(id_cell!=quadcell)
      {
        id_new_cell = dst->InsertNextCell(type_cell, N_pts, pts);
        copyCellData(src, id_cell, dst, id_new_cell);
      }
      else
      {
        vtkIdType triangle1[3];
        vtkIdType triangle2[3];
        triangle1[0]=pts[1];
        triangle1[1]=pts[3];
        triangle1[2]=pts[0];
        triangle2[0]=pts[3];
        triangle2[1]=pts[1];
        triangle2[2]=pts[2];
        id_new_cell1 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle1);
        copyCellData(src, id_cell, dst, id_new_cell1);
        id_new_cell2 = dst->InsertNextCell(VTK_TRIANGLE, 3, triangle2);
        copyCellData(src, id_cell, dst, id_new_cell2);
        stencil_t S;
        S.id_cell1=id_new_cell1;
        S.id_cell2=id_new_cell2;
        S.p[0]=pts[0];
        S.p[1]=pts[1];
        S.p[2]=pts[2];
        S.p[3]=pts[3];
        S.valid=true;
        SwapCell(dst,S);
      }
    };
    cout_grid(cout,dst,true,true,true,true);
    makeCopy(dst, src);
  }//end of if quad
}

void DeletePickedCell::operate()
{
  vtkIdType cellId = GuiMainWindow::pointer()->getPickedCell();
  cout<<"You picked "<<cellId<<endl;
  foobar(grid,cellId);
};
