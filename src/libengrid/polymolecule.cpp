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

#include "polymolecule.h"

PolyMolecule::PolyMolecule(PolyMesh *poly_mesh, const QList<int> &faces)
{
  m_PolyMesh = poly_mesh;
  QSet<int> nodes;
  QHash<int,int> node_map;
  foreach (int face, faces) {
    for (int i_node = 0; i_node < m_PolyMesh->numNodes(face); ++i_node) {
      nodes.insert(m_PolyMesh->nodeIndex(face, i_node));
    }
  }
  {
    int N = 0;
    foreach (int node, nodes) {
      m_Nodes.append(node);
      node_map[N] = node;
      ++N;
    }
  }
  m_Faces.resize(faces.size());
  for (int i_face = 0; i_face < faces.size(); ++i_face) {
    for (int i_node = 0; i_node < m_PolyMesh->numNodes(faces[i_face]); ++i_node) {
      m_Faces[i_face].append(node_map[m_PolyMesh->nodeIndex(faces[i_face], i_node)]);
    }
    m_Faces[i_face].append(m_Faces[i_face].first());
  }
}

vec3_t PolyMolecule::getXFace(int face)
{

}
