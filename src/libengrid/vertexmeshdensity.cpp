// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "vertexmeshdensity.h"
#include "egvtkobject.h"
#include <iostream>
using namespace std;

VertexMeshDensity::VertexMeshDensity()
{
  type = EG_SIMPLE_VERTEX;
  CurrentNode = 0;
}

//this=user defined
//VMD=VMD of current node
//This operator is NOT SYMMETRICAL. But it can be used with indexOf.
bool VertexMeshDensity::operator==(const VertexMeshDensity & VMD) const
{
  if(this->nodeset.size()>0) //node ID mode
  {
    if(this->nodeset.contains(VMD.CurrentNode) ) return(true);//node ID given, we're done
    else return(false);
  }
  else //node properties mode
  {
    if(this->type!=VMD.type && this->type!=-1) return(false);
    
    QMapIterator<int, int> i(this->BCmap);
    while (i.hasNext()) {
      i.next();
      int index=i.key();
      if((this->BCmap)[index]!=VMD.BCmap[index] && (this->BCmap)[index]!=1) return(false);
    }
    
    return(true);
  }
}

int VertexMeshDensity::findSmallestVMD( QVector <VertexMeshDensity> vector)
{
  int ret = -1;
  double Lmin = 0;
  bool first = true;
  for(int i=0;i<vector.size();i++){
    if(vector[i]==*this) {
      if(first) {
        ret = i;
        Lmin = vector[i].density;
        first = false;
      }
      else if(vector[i].density < Lmin) {
        ret = i;
        Lmin = vector[i].density;
      }
    }
  }
  return(ret);
}

//converts string to nodeset
void VertexMeshDensity::setNodes(QString str)
{
  nodeset.clear();//empty by default
//   cout<<"str.size="<<str.size()<<endl;
  if(str.size()>0)
  {
    QStringList L = str.split(",");
    foreach(QString elt,L) nodeset.insert(elt.toInt());
  }
}

ostream& operator<<(ostream &out, VertexMeshDensity A)
{
  out<<" BCmap="<<A.BCmap;
  out<<" type="<<(int)A.type;
  out<<" nodeset="<<A.nodeset;
  out<<" density="<<A.density;
  return(out);
}

ostream& operator<<(ostream &out, QVector<VertexMeshDensity> VMDvector)
{
  int N=VMDvector.size();
  out<<"Size="<<N;
  if(N>0) out<<endl;
  for(int i=0;i<N;i++)
  {
    out<<VMDvector[i];
    if(i!=N-1) out<<endl;
  }
  return(out);
}
