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
  if(this->nodeset.contains(VMD.CurrentNode) ) return(true);
//   if(this->type==VMD.type && this->BClist==VMD.BClist) return(true);
  else return(false);
}

void VertexMeshDensity::SetNodes(QString str)
{
  QStringList L = str.split(",");
  nodeset.clear();
  foreach(QString elt,L) nodeset.insert(elt.toInt());
}

// QVector < QVector <int> > BClist;
// QVector <char> type;
// QVector <double> density;

ostream& operator<<(ostream &out, VertexMeshDensity A)
{
  out<<" BClist="<<A.BClist;
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
