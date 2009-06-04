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
#include "showinfo.h"

#include "insertpoints.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <vtkCharArray.h>

ShowInfo::ShowInfo(bool b, vtkIdType P, vtkIdType C) : SurfaceOperation()
{
  CellInfo = b;
  PickedPoint = P;
  PickedCell = C;
}

void ShowInfo::operate()
{
  int N_cells=grid->GetNumberOfCells();
  int N_points=grid->GetNumberOfPoints();
  if(CellInfo)
  {
    if(PickedCell>=0 && PickedCell<N_cells)
    {
      cout<<"=== INFO ON CELL "<<PickedCell<<" ==="<<endl;
      QVector<vtkIdType> absolute_c2c;
      foreach(int i, c2c[_cells[PickedCell]]){
        if(i!=-1) absolute_c2c.push_back(cells[i]);
      }
      
      cout<<"absolute_c2c(PickedCell)="<<absolute_c2c<<endl;
      vtkIdType *pts, N_pts;
      grid->GetCellPoints(PickedCell, N_pts, pts);
      cout<<"pts=";
      for(int i=0;i<N_pts;i++) cout<<pts[i]<<" ";
      cout<<endl;
      cout<<"coords:"<<endl;
      for(int i=0;i<N_pts;i++) {
        vec3_t X;
        grid->GetPoint(pts[i],X.data());
        cout<<"pts["<<i<<"]="<<X<<endl;
      }
      cout<<"area="<<cellVA(grid,PickedCell)<<endl;
      cout<<"Q_L("<<PickedCell<<")="<<Q_L(PickedCell)<<endl;
      cout<<"====================================="<<endl;
    }
    else
    {
      cout<<"Invalid cell"<<endl;
    }
  }
  else
  {
    if(PickedPoint>=0 && PickedPoint<N_points)
    {
      cout<<"=== INFO ON POINT "<<PickedPoint<<" ==="<<endl;
      
      QSet<vtkIdType> absolute_n2n;
      foreach(int i, n2n[_nodes[PickedPoint]]){
        if(i!=-1) absolute_n2n.insert(nodes[i]);
      }
      cout<<"absolute_n2n(PickedPoint)="<<absolute_n2n<<endl;
      
      QSet<vtkIdType> absolute_n2c;
      foreach(int i, n2c[_nodes[PickedPoint]]){
        if(i!=-1) absolute_n2c.insert(cells[i]);
      }
      cout<<"absolute_n2c(PickedPoint)="<<absolute_n2c<<endl;
      
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      cout<<"node_type="<<VertexType2Str(node_type->GetValue(PickedPoint))<<endl;
      vec3_t X;
      grid->GetPoint(PickedPoint,X.data());
      cout<<"X="<<X<<endl;
      cout<<"G_k("<<PickedPoint<<")="<<G_k(PickedPoint)<<endl;
      cout<<"====================================="<<endl;
    }
    else
    {
      cout<<"Invalid point"<<endl;
    }
  }
}
