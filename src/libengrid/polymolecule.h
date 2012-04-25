//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                      +
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

#ifndef POLYMOLECULE_H
#define POLYMOLECULE_H

#include "egvtkobject.h"
#include "polymesh.h"

#include <QList>

class PolyMolecule : public EgVtkObject
{

private: // data types

  struct edge_t {
    int node1;
    int node2;
    int face1;
    int face2;
  };

private: // attributes

  PolyMesh*             m_PolyMesh;
  QList<int>            m_Nodes;
  QVector<QList<int> >  m_Faces;
  QList<QList<edge_t> > m_N2N;

public:

  PolyMolecule(PolyMesh *poly_mesh, const QList<int> &faces);

  vec3_t getXNode(int node) { return m_PolyMesh->nodeVector(m_Nodes[node]); }
  vec3_t getXFace(int face);

};

#endif // POLYMOLECULE_H
