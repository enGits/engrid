//
// C++ Implementation: showinfo
//
// Description: 
//
//
// Author: Mike Taverne <mtaverne@engits.com>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "showinfo.h"

#include "geometrytools.h"
using namespace GeometryTools;

#include <vtkCharArray.h>

ShowInfo::ShowInfo(bool b, vtkIdType P, vtkIdType C) : Operation()
{
  CellInfo=b;
  PickedPoint=P;
  PickedCell=C;
}


// ShowInfo::~ShowInfo()
// {
// }


void ShowInfo::operate()
{
  int N_cells=grid->GetNumberOfCells();
  int N_points=grid->GetNumberOfPoints();
  if(CellInfo)
  {
    if(PickedCell>=0 && PickedCell<N_cells)
    {
      cout<<"=== INFO ON CELL "<<PickedCell<<" ==="<<endl;
      cout<<"c2c_func(PickedCell)="<<c2c_func(PickedCell)<<endl;
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
      cout<<"n2n_func(PickedPoint)="<<n2n_func(PickedPoint)<<endl;
      cout<<"n2c_func(PickedPoint)="<<n2c_func(PickedPoint)<<endl;
      EG_VTKDCN(vtkCharArray, node_type, grid, "node_type");//node type
      cout<<"node_type="<<VertexType2Str(node_type->GetValue(PickedPoint))<<endl;
      vec3_t X;
      grid->GetPoint(PickedPoint,X.data());
      cout<<"X="<<X<<endl;
      cout<<"====================================="<<endl;
    }
    else
    {
      cout<<"Invalid point"<<endl;
    }
  }
}
