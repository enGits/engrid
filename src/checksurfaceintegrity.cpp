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
// #include "egvtkobject.h"

CheckSurfaceIntegrity::CheckSurfaceIntegrity()
 : Operation()
{
}

void CheckSurfaceIntegrity::operate()
{
  cout<<"this->isWaterTight()="<<this->isWaterTight()<<endl;
  cout<<"this->Nmin="<<this->Nmin<<endl;
  cout<<"this->Nmax="<<this->Nmax<<endl;
  cout<<"this->BadCells="<<this->BadCells<<endl;
}

bool CheckSurfaceIntegrity::isWaterTight()
{
  setAllSurfaceCells();
  BadCells.clear();
  
  bool first = true;
  foreach(vtkIdType node1,nodes) {
    foreach(vtkIdType node2,n2n_func(node1)) {
      QSet <vtkIdType> edge_cells;
      int N = getEdgeCells(node1,node2,edge_cells);
      if(first) {
        first = false;
        Nmin = N;
        Nmax = N;
      }
      else {
        Nmin = min(Nmin,N);
        Nmax = max(Nmax,N);
      }
      if(edge_cells.size()!=2) BadCells.unite(edge_cells);
    }
  }
  if( Nmin==2 && Nmax==2 ) return(true);
  else return(false);
}
