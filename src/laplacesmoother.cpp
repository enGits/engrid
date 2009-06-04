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
#include "laplacesmoother.h"
#include <vtkCellLocator.h>
#include <vtkCharArray.h>
#include <vtkGenericCell.h>
#include "guimainwindow.h"

using namespace GeometryTools;

LaplaceSmoother::LaplaceSmoother()
: SurfaceOperation()
{
   DebugLevel=0;
   setQuickSave(true);
}

void LaplaceSmoother::operate()
{
  cout<<"=== LaplaceSmoother START ==="<<endl;
  
  QVector<vtkIdType> AllCells;
  getAllSurfaceCells(AllCells, grid);
  QVector<vtkIdType> SelectedCells;
  getSurfaceCells(m_bcs, SelectedCells, grid);
  
  cout<<"setCells START"<<endl;
  setCells(AllCells);
  cout<<"setCells END"<<endl;
  
  QSet <vtkIdType> SelectedNodes;
  getSurfaceNodes(m_bcs,SelectedNodes,grid);
  
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");
  int moved_points=0;
  
  for(int i_iter=0;i_iter<NumberOfIterations;i_iter++)
  {
    foreach(vtkIdType id_nodeG,SelectedNodes)
    {
      if(node_type->GetValue(id_nodeG)==VTK_SIMPLE_VERTEX)
      {
        vec3_t G(0,0,0);
        QVector <vtkIdType> PSP = getPotentialSnapPoints(id_nodeG);
        foreach(vtkIdType id_nodeM, PSP)
        {
          vec3_t M;
          grid->GetPoint(id_nodeM, M.data());
          G+=M;
        }
        G=(1./PSP.size())*G;
        
        vec3_t P;
        P=project(G);
        
        //check that no cell gets flipped!
        while(FlippedCells(id_nodeG,P))
        {
          vec3_t x0_old;
          grid->GetPoint(id_nodeG, x0_old.data());
          P=x0_old+0.5*(P-x0_old);
        };
        
        grid->GetPoints()->SetPoint(id_nodeG, P.data());
        
        moved_points++;
      }
    }
  }
  
  cout<<"=== LaplaceSmoother END ==="<<endl;
}

bool LaplaceSmoother::FlippedCells(vtkIdType id_nodeG, vec3_t P)
{
  vec3_t x0_old, x0_new;
  grid->GetPoint(id_nodeG, x0_old.data());
  x0_new=P;
  
  foreach(int i_cell, n2c[_nodes[id_nodeG]])
  {
    vtkIdType id_cell = cells[i_cell];
    vtkIdType N_pts, *pts;
    grid->GetCellPoints(id_cell, N_pts, pts);
    int i;
    for(i=0;i<N_pts;i++)
    {
      if(pts[i]==id_nodeG) break;
    }
    vec3_t x2, x3;
    grid->GetPoint(pts[(i+1)%N_pts], x2.data());
    grid->GetPoint(pts[(i+2)%N_pts], x3.data());
    vec3_t v2_old=x2-x0_old;
    vec3_t v3_old=x3-x0_old;
    
    //top point
    vec3_t S=v2_old.cross(v3_old);
    double V_old=tetraVol(x0_old, S, x2, x3, true);
    double V_new=tetraVol(x0_new, S, x2, x3, true);
    double prod=V_old*V_new;
    if( prod<0 ) {
      return(true);
    }
  }
  return(false);
}
