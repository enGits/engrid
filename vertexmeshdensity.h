#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>
#include <QVector>

class VertexMeshDensity{
public:
  VertexMeshDensity();
  ~VertexMeshDensity();

public:
  QVector < QVector <int> > BClist;
  QVector <char> type;
  QVector <double> density;
};

ostream& operator<<(ostream &out, VertexMeshDensity A);

#endif
