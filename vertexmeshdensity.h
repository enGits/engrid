#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>

class VertexMeshDensity{
public:
  VertexMeshDensity();
  ~VertexMeshDensity();

public:
  QSet <vtkIdType> NeighbourSet;
  char type;
};

#endif
