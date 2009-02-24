#include "vertexmeshdensity.h"

VertexMeshDensity::VertexMeshDensity()
{
}


VertexMeshDensity::~VertexMeshDensity()
{
}

// QVector < QVector <int> > BClist;
// QVector <char> type;
// QVector <double> density;

ostream& operator<<(ostream &out, VertexMeshDensity A)
{
  int N1=A.type.size();
  out<<"Size="<<N1<<endl;
  for(int i=0;i<N1;i++){
    int N2=A.BClist[i].size();
    out<<"BClist=[";
    for(int j=0;j<N2;j++){
      out<<A.BClist[i][j];
      if(j!=N2-1) out<<",";
    }
    out<<"]";
    out<<" type="<<(int)A.type[i]<<" density="<<A.density[i];
    if(i!=N1-1) out<<endl;
  }
  return(out);
}
