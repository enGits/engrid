#ifndef VERTEXMESHDENSITY_H
#define VERTEXMESHDENSITY_H

#include <vtkIdList.h>
#include <QSet>
#include <QVector>
#include <QString>
#include "egvtkobject.h"

class VertexMeshDensity{
public:
  VertexMeshDensity();
//   ~VertexMeshDensity();

public:
  QVector <int> BClist;
  QVector <Qt::CheckState> BClist_value;
  char type;
  QSet <vtkIdType> nodeset;
  double density;
  vtkIdType CurrentNode;
  bool operator==(const VertexMeshDensity & VMD) const;
  void SetNodes(QString str);
};

ostream& operator<<(ostream &out, VertexMeshDensity A);

ostream& operator<<(ostream &out, QVector<VertexMeshDensity> VMDvector);

#endif
