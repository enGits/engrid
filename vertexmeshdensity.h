#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>
#include <QVector>

class VertexMeshDensity{
public:
  VertexMeshDensity();
//   ~VertexMeshDensity();

public:
  QVector <int> BClist;
  char type;
  double density;
  bool operator==(const VertexMeshDensity & VMD) const;
};

ostream& operator<<(ostream &out, VertexMeshDensity A);

ostream& operator<<(ostream &out, QVector<VertexMeshDensity> VMDvector);

#endif
