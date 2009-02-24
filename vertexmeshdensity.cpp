#include "vertexmeshdensity.h"

VertexMeshDensity::VertexMeshDensity()
{
  NeighbourSet=QSet <vtkIdType> ();
  type=0;
}


VertexMeshDensity::~VertexMeshDensity()
{
}
