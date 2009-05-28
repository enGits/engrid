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
#include "checksurfaceintegrity.h"

CheckSurfaceIntegrity::CheckSurfaceIntegrity()
 : Operation()
{
}

void CheckSurfaceIntegrity::operate()
{
  cout<<"this->isWaterTight()="<<this->isWaterTight()<<endl;
}

bool CheckSurfaceIntegrity::isWaterTight()
{
  setAllSurfaceCells();
  bool first = true;
  foreach(vtkIdType node1,nodes) {
    foreach(vtkIdType node2,n2n_func(node1)) {
      QVector <vtkIdType> edge_cells = getEdgeCells(node1,node2);
      if(first) {
        first = false;
        Nmin = edge_cells.size();
        Nmax = edge_cells.size();
      }
      else {
        Nmin = min(Nmin,edge_cells.size());
        Nmax = max(Nmax,edge_cells.size());
      }
    }
  }
  if( Nmin==2 && Nmax==2 ) return(true);
  else return(false);
}
