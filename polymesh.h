//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008 Oliver Gloth                                          +
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
#ifndef polymesh_H
#define polymesh_H

class PolyMesh;

#include "egvtkobject.h"

class PolyMesh : public EgVtkObject
{
  
protected: // data types
  
  enum idx_t { node, edge, face, cell };
  
  struct node_t {
    PolyMesh *poly; // 4
    idx_t type;     // 4
    vtkIdType idx;  // 4
    int subidx;     // 4
    //double w;       // 8
    //vec3_t x;       // 24
    //vec3_t n;       // 24
    float x,y,z;
    node_t(PolyMesh *p, idx_t t, int i, int si) { poly = p; type = t; idx = i; subidx = si; };
    node_t() { poly = NULL, type = node, idx = 0; subidx = 0; };
    bool operator==(const node_t N) const;
  };
  
  friend uint qHash(node_t N);
  friend ostream& operator<<(ostream& s, node_t N);
  
  struct face_t {
    QVector<node_t> node;
    int owner, neighbour, bc;
    node_t operator[](int i) { while(i<0) i+=node.size(); while(i>=node.size()) i-=node.size(); return node[i]; };
    void checkOrientation();
    bool operator<(const face_t &F) const;
    vec3_t normal();
  };
  
  friend ostream& operator<<(ostream& s, face_t F);

protected: // attributes
  
  bool dbg;
  vtkUnstructuredGrid *grid;
  QList<face_t> face_list;
  QVector<face_t> faces;
  QVector<int> cell2pc;
  QVector<int> node2pc;
  int id_newpc;
  QVector<vtkIdType> cells;
  QVector<vec3_t> nodes;
  QVector<QVector<int> > c2c;
  QList<int> bcs;
  QVector<double> weight;
  bool dual;
  
protected: // methods
  
  int pcIdxNode(vtkIdType id_node);
  int pcIdxCell(vtkIdType id_node);
  void pass1Tetras();
  void pass1Prisms();
  void pass1Hexas();
  void pass1();
  void computeNodes();
  void createNodes();
  void sortFaces();
  face_t combineFaces(QList<face_t> faces);
  double faceW(double w) { if (w > 1.01) return 0.5*w; else return w; };
  
  void pass2();
  void pass3();
  
public: // methods
  
  PolyMesh(vtkUnstructuredGrid *a_grid, bool dual_mesh = true);
  int totalNumNodes() const { return nodes.size(); };
  vec3_t nodeVector(int i) const { return nodes[i]; };
  int numNodes(int i) const { return faces[i].node.size(); };
  int nodeIndex(int i, int j) const { return faces[i].node[j].idx; };
  int numFaces() const { return faces.size(); };
  int numCells() const { return id_newpc; };
  int owner(int i) const { return faces[i].owner; };
  int neighbour(int i) const { return faces[i].neighbour; };
  int boundaryCode(int i) const { return faces[i].bc; };
  int numBCs() const { return bcs.size()-1; };
  void getFace(vtkIdType idx, int subidx, QList<vtkIdType> &nodes);
  void getEdge(vtkIdType idx, int subidx, QList<vtkIdType> &nodes);
  
  
};

uint qHash(PolyMesh::node_t N);

#endif