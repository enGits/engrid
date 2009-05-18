#include "vertexmeshdensity.h"
#include "egvtkobject.h"
#include <iostream>
using namespace std;

VertexMeshDensity::VertexMeshDensity()
{
}


// VertexMeshDensity::~VertexMeshDensity()
// {
// }

//this=user defined
//VMD=VMD of current node
//This operator is NOT SYMMETRICAL. But it can be used with indexOf.
bool VertexMeshDensity::operator==(const VertexMeshDensity & VMD) const
{
/*  cout<<"this->nodeset="<<this->nodeset<<endl;
  cout<<"VMD.nodeset="<<VMD.nodeset<<endl;*/
  
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
//     cout << i.key() << ": " << i.value() << endl;
      int index=i.key();
      if((this->BCmap)[index]!=VMD.BCmap[index] && (this->BCmap)[index]!=1) return(false);
    }
    
    return(true);
  }
}

//converts string to nodeset
void VertexMeshDensity::SetNodes(QString str)
{
  nodeset.clear();//empty by default
  cout<<"str.size="<<str.size()<<endl;
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
