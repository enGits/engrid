//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2013 enGits GmbH                                      +
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
  QList<int>            m_FacesInPMesh;
  QList<int>            m_Nodes;
  QVector<QList<int> >  m_Faces;
  QVector<vec3_t>       m_FaceNormals;
  QVector<vec3_t>       m_NodeNormals;
  QVector<QSet<int> >   m_N2N;
  QVector<QSet<int> >   m_F2F;
  QVector<QSet<int> >   m_N2BadFaces;
  PolyMolecule*         m_SplitCell1;
  PolyMolecule*         m_SplitCell2;

  int                   m_PCell;
  vec3_t                m_CentreOfGravity;
  double                m_MaxPyramidVolume;
  double                m_MinPyramidVolume;
  bool                  m_AllPositive;

private: // methods

  void computeCentreOfGravity();
  void buildNode2Node();
  void buildFace2Face();
  void computeNormals();
  void smooth(bool delaunay = true, bool write = false);
  void split(bool write = false);
  void updateFace(int face, int new_cell_index);
  void centreSplit();


  template <class C> void init(PolyMesh *poly_mesh, const C &faces);
  template <class C> void setSubMolecules(const C &face_indices1);

  edge_t getEdge(int node1, int node2);
  QList<edge_t> findConcaveEdges();

public:

  PolyMolecule();
  PolyMolecule(PolyMesh *poly_mesh, int i_pcell);

  vec3_t getXNode(int node) { return m_PolyMesh->nodeVector(m_Nodes[node]); }
  vec3_t getXFace(int face);
  void   writeVtkFile(QString file_name);
  void   createPolyData(vtkPolyData *poly_data);
  double minPyramidVolume() { return m_MinPyramidVolume; }
  double maxPyramidVolume() { return m_MaxPyramidVolume; }
  void   fix(bool write = false);
  bool   allPositive() { return m_AllPositive; }
  void   updatePMesh();

};

template <class C>
void PolyMolecule::init(PolyMesh *poly_mesh, const C &faces)
{
  m_PolyMesh = poly_mesh;
  QSet<int> nodes;
  QHash<int,int> node_map;
  for (int i = 0; i < faces.size(); ++i) {
    for (int j = 0; j < faces[i].size(); ++j) {
      nodes.insert(faces[i][j]);
    }
  }
  {
    int N = 0;
    foreach (int node, nodes) {
      m_Nodes.append(node);
      node_map[node] = N;
      ++N;
    }
  }
  m_Faces.resize(faces.size());
  for (int i = 0; i < faces.size(); ++i) {
    for (int j = 0; j < faces[i].size(); ++j) {
      m_Faces[i].append(node_map[faces[i][j]]);
    }
  }
  if (m_Nodes.size() == 0) {
    EG_BUG;
  }
  if (m_Faces.size() == 0) {
    EG_BUG;
  }
  computeCentreOfGravity();
  buildNode2Node();
  computeNormals();
}

template <class C>
void PolyMolecule::setSubMolecules(const C &face_indices1)
{
  QList<int> face_indices2;
  for (int i = 0; i < m_Faces.size(); ++i) {
    if (!face_indices1.contains(i)) {
      face_indices2.append(i);
    }
  }
  QList<QList<int> > faces1;
  foreach (int i, face_indices1) {
    QList<int> face;
    foreach (int node, m_Faces[i]) {
      face.append(node);
    }
    faces1.append(face);
  }
  delete m_SplitCell1;
  m_SplitCell1 = new PolyMolecule();
  m_SplitCell1->init(m_PolyMesh, faces1);
  QList<QList<int> > faces2;
  foreach (int i, face_indices2) {
    QList<int> face;
    foreach (int node, m_Faces[i]) {
      face.append(node);
    }
    faces2.append(face);
  }
  delete m_SplitCell2;
  m_SplitCell2 = new PolyMolecule();
  m_SplitCell2->init(m_PolyMesh, faces2);
}

#endif // POLYMOLECULE_H
