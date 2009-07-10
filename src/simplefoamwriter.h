//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008,2009 Oliver Gloth                                     +
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
#ifndef simplefoamwriter_H
#define simplefoamwriter_H

class SimpleFoamWriter;

#include "iooperation.h"

/**
 * Writer for OpenFOAM grids
 */
class SimpleFoamWriter : public IOOperation
{

protected: // data types
  
  struct face_t {
    QVector<vtkIdType> node;
    vtkIdType owner, neighbour;
    int bc;
    int operator[](int i) { while(i<0) i+=node.size(); while(i>=node.size()) i-=node.size(); return node[i]; }
    bool operator<(const face_t &F) const;
    vec3_t normal(vtkUnstructuredGrid *grid);
    face_t() {}
    face_t(int N, int o, int n, int b=0) { node.resize(N); owner = o; neighbour = n; bc = b; }
  };

  struct patch_t {
    int bc, startFace, nFaces;
    bool operator<(const patch_t &P) const { return startFace < P.startFace; }
  };
  
protected: // attributes
  
  QString         m_Path;
  QVector<face_t> m_Faces;
  QList<face_t>   m_LFaces;
  vtkIntArray*    m_BC;
  int             m_NumCells;
  QVector<int>    m_Eg2Of;
  
protected: // methods
  
  vtkIdType getNeigh(int i_cells, int i_neigh);
  void addFace(face_t F);
  void createFaces();
  void writePoints();
  void writeFaces();
  void writeOwner();
  void writeNeighbour();
  void writeBoundary();
  
  virtual void operate();
  
public: // methods
  
  SimpleFoamWriter();
  virtual ~SimpleFoamWriter() {}
  
};

#endif
