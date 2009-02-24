#include "vertexmeshdensity.h"

VertexMeshDensity::VertexMeshDensity()
{
}


// VertexMeshDensity::~VertexMeshDensity()
// {
// }

bool VertexMeshDensity::operator==(const VertexMeshDensity & VMD) const
{
  if(this->type==VMD.type && this->BClist==VMD.BClist) return(true);
  else return(false);
}

// QVector < QVector <int> > BClist;
// QVector <char> type;
// QVector <double> density;

ostream& operator<<(ostream &out, VertexMeshDensity A)
{
  int N=A.BClist.size();
  out<<"BClist=[";
  for(int i=0;i<N;i++){
    out<<A.BClist[i];
    if(i!=N-1) out<<",";
  }
  out<<"]";
  out<<" type="<<(int)A.type<<" density="<<A.density;
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
