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

CheckSurfaceIntegrity::CheckSurfaceIntegrity() : SurfaceOperation()
{
  EG_TYPENAME;
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

  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2l_t  n2n   = getPartN2N();
  
  bool first = true;
  foreach(vtkIdType id_node1, nodes) {
    foreach(int i_node2, n2n[_nodes[id_node1]]) {
      vtkIdType id_node2 = nodes[i_node2];
      QSet <vtkIdType> edge_cells;
      int N = getEdgeCells(id_node1, id_node2,edge_cells);
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
