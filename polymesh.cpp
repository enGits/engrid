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
#include "polymesh.h"

uint qHash(PolyMesh::node_t N) 
{ 
  uint h = 0;
  if (N.type == PolyMesh::node) {
    h = N.idx;
  } else {
    if (N.type == PolyMesh::cell) {
      vtkIdType N_pts, *pts;
      N.poly->grid->GetCellPoints(N.idx,N_pts,pts);
      for (int i = 0; i < N_pts; ++i) h += pts[i];
    } else {
      QList<vtkIdType> nds;
      if (N.type == PolyMesh::edge) N.poly->getEdge(N.idx,N.subidx,nds);
      if (N.type == PolyMesh::face) N.poly->getFace(N.idx,N.subidx,nds);
      foreach (vtkIdType p, nds) h += p;
    };
  };
  return h;
};


bool PolyMesh::node_t::operator==(const node_t N) const
{
  QList<vtkIdType> nds;
  if (type == edge) {
    poly->getEdge(idx,subidx,nds);
  };
  if (type == face) {
    poly->getFace(idx,subidx,nds);
  };
  bool match = ((type==N.type) && (idx==N.idx) && (subidx==N.subidx)); 
  if (!match && (type == N.type) && (idx != N.idx)) {
    if ((type == edge) || (type == face)) {
      match = true;
      QList<vtkIdType> nds1, nds2;
      if (type == edge) {
        poly->getEdge(idx,subidx,nds1);
        poly->getEdge(N.idx,N.subidx,nds2);
      };
      if (type == face) {
        poly->getFace(idx,subidx,nds1);
        poly->getFace(N.idx,N.subidx,nds2);
      };
      foreach (vtkIdType p, nds1) {
        if (!nds2.contains(p)) {
          match = false;
          break;
        };
      };
    };
  };
  return match;
};

ostream& operator<<(ostream &s, PolyMesh::node_t N)
{
  if (N.type == PolyMesh::node) s << "node: ";
  if (N.type == PolyMesh::edge) s << "edge: ";
  if (N.type == PolyMesh::face) s << "face: ";
  if (N.type == PolyMesh::cell) s << "cell: ";
  s << N.idx << ',' << N.subidx;
  if (N.type != PolyMesh::node) {
    s << "(";
    if (N.type == PolyMesh::cell) {
      vtkIdType N_pts, *pts;
      N.poly->grid->GetCellPoints(N.idx,N_pts,pts);
      for (int i = 0; i < N_pts-1; ++i) s << pts[i] << ',';
      s << pts[N_pts-1];
    } else {
      QList<vtkIdType> nds;
      if (N.type == PolyMesh::edge) {
        N.poly->getEdge(N.idx,N.subidx,nds);
      } else {
        N.poly->getFace(N.idx,N.subidx,nds);
      };
      QList<vtkIdType>::iterator i = nds.begin();
      s << *i;
      ++i;
      while (i != nds.end()) {
        s << ',' << *i;
        ++i;
      };
    }; 
    s << ")";
  };
  return s;
};

ostream& operator<<(ostream &s, PolyMesh::face_t F)
{
  s << "bc: " << F.bc << " owner: " << F.owner << " neighbour: " << F.neighbour;
  for (int i = 0; i < F.node.size(); ++i) {
    s << "\n  " << F.node[i];
  };
  return s;
};

void PolyMesh::face_t::checkOrientation()
{
  if (owner > neighbour) {
    QVector<node_t> node_copy(node.size());
    qCopy(node.begin(),node.end(),node_copy.begin());
    for (int i = 0; i < node.size(); ++i) {
      node[i] = node_copy[node.size()-i-1];
    };
    swap(owner,neighbour);
  };
};

bool PolyMesh::face_t::operator<(const face_t &F) const
{
  bool less = false;
  if (bc < F.bc) {
    less = true;
  } else if (bc == F.bc) {
    if (owner < F.owner) {
      less = true;
    } else if (owner == F.owner) {
      if (neighbour < F.neighbour) {
        less = true;
      };
    };
  };
  
  return less;
};

vec3_t PolyMesh::face_t::normal()
{
  vec3_t x(0,0,0);
  if (node.size() < 3) EG_BUG;
  for (int i = 0; i < node.size(); ++i) {
    x += vec3_t(node[i].x,node[i].y,node[i].z);
  };
  x *= 1.0/node.size();
  vec3_t n(0,0,0);
  for (int i = 0; i < node.size()-1; ++i) {
    vec3_t a = vec3_t(node[i].x,node[i].y,node[i].z)-x; 
    vec3_t b = vec3_t(node[i+1].x,node[i+1].y,node[i+1].z)-x; 
    n += 0.5*(a.cross(b));
  };
  vec3_t a = vec3_t(node[node.size()-1].x,node[node.size()-1].y,node[node.size()-1].z) - x; 
  vec3_t b = vec3_t(node[0].x,node[0].y,node[0].z) - x; 
  n += 0.5*(a.cross(b));
  return n;
};

void PolyMesh::getFace(vtkIdType idx, int subidx, QList<vtkIdType> &nodes)
{
  nodes.clear();
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(idx, N_pts, pts);
  vtkIdType type_cell = grid->GetCellType(idx);
  if (type_cell == VTK_TETRA) {
    if (subidx == 0) {
      nodes.append(pts[2]);
      nodes.append(pts[1]);
      nodes.append(pts[0]);
    };
    if (subidx == 1) {
      nodes.append(pts[1]);
      nodes.append(pts[3]);
      nodes.append(pts[0]);
    };
    if (subidx == 2) {
      nodes.append(pts[3]);
      nodes.append(pts[2]);
      nodes.append(pts[0]);
    };
    if (subidx == 3) {
      nodes.append(pts[2]);
      nodes.append(pts[3]);
      nodes.append(pts[1]);
    };
  };
  if (type_cell == VTK_WEDGE) {
    if (subidx == 0) {
      nodes.append(pts[0]);
      nodes.append(pts[1]);
      nodes.append(pts[2]);
    };
    if (subidx == 1) {
      nodes.append(pts[3]);
      nodes.append(pts[5]);
      nodes.append(pts[4]);
    };
    if (subidx == 2) {
      nodes.append(pts[3]);
      nodes.append(pts[4]);
      nodes.append(pts[1]);
      nodes.append(pts[0]);
    };
    if (subidx == 3) {
      nodes.append(pts[1]);
      nodes.append(pts[4]);
      nodes.append(pts[2]);
      nodes.append(pts[5]);
    };
    if (subidx == 4) {
      nodes.append(pts[0]);
      nodes.append(pts[2]);
      nodes.append(pts[5]);
      nodes.append(pts[3]);
    };
  };
  if (type_cell == VTK_HEXAHEDRON) {
    if (subidx == 0) {
      nodes.append(pts[0]);
      nodes.append(pts[3]);
      nodes.append(pts[2]);
      nodes.append(pts[1]);
    };
    if (subidx == 1) {
      nodes.append(pts[4]);
      nodes.append(pts[5]);
      nodes.append(pts[6]);
      nodes.append(pts[7]);
    };
    if (subidx == 2) {
      nodes.append(pts[0]);
      nodes.append(pts[1]);
      nodes.append(pts[5]);
      nodes.append(pts[4]);
    };
    if (subidx == 3) {
      nodes.append(pts[3]);
      nodes.append(pts[7]);
      nodes.append(pts[6]);
      nodes.append(pts[2]);
    };
    if (subidx == 4) {
      nodes.append(pts[0]);
      nodes.append(pts[4]);
      nodes.append(pts[7]);
      nodes.append(pts[3]);
    };
    if (subidx == 5) {
      nodes.append(pts[1]);
      nodes.append(pts[2]);
      nodes.append(pts[6]);
      nodes.append(pts[5]);
    };
  };
};

void PolyMesh::getEdge(vtkIdType idx, int subidx, QList<vtkIdType> &nodes)
{
  nodes.clear();
  vtkIdType N_pts, *pts;
  grid->GetCellPoints(idx, N_pts, pts);
  vtkIdType type_cell = grid->GetCellType(idx);
  if (type_cell == VTK_TETRA) {
    if (subidx == 0) {
      nodes.append(pts[0]);
      nodes.append(pts[1]);
    };
    if (subidx == 1) {
      nodes.append(pts[0]);
      nodes.append(pts[2]);
    };
    if (subidx == 2) {
      nodes.append(pts[0]);
      nodes.append(pts[3]);
    };
    if (subidx == 3) {
      nodes.append(pts[1]);
      nodes.append(pts[2]);
    };
    if (subidx == 4) {
      nodes.append(pts[1]);
      nodes.append(pts[3]);
    };
    if (subidx == 5) {
      nodes.append(pts[2]);
      nodes.append(pts[3]);
    };
  };
  if (type_cell == VTK_WEDGE) {
    if (subidx == 0) {
      nodes.append(pts[0]);
      nodes.append(pts[1]);
    };
    if (subidx == 1) {
      nodes.append(pts[0]);
      nodes.append(pts[2]);
    };
    if (subidx == 2) {
      nodes.append(pts[0]);
      nodes.append(pts[3]);
    };
    if (subidx == 3) {
      nodes.append(pts[1]);
      nodes.append(pts[2]);
    };
    if (subidx == 4) {
      nodes.append(pts[1]);
      nodes.append(pts[4]);
    };
    if (subidx == 5) {
      nodes.append(pts[2]);
      nodes.append(pts[5]);
    };
    if (subidx == 6) {
      nodes.append(pts[3]);
      nodes.append(pts[4]);
    };
    if (subidx == 7) {
      nodes.append(pts[3]);
      nodes.append(pts[5]);
    };
    if (subidx == 8) {
      nodes.append(pts[4]);
      nodes.append(pts[5]);
    };
  };
  if (type_cell == VTK_HEXAHEDRON) {
    if (subidx == 0) {
      nodes.append(pts[0]);
      nodes.append(pts[1]);
    };
    if (subidx == 1) {
      nodes.append(pts[0]);
      nodes.append(pts[3]);
    };
    if (subidx == 2) {
      nodes.append(pts[0]);
      nodes.append(pts[4]);
    };
    if (subidx == 3) {
      nodes.append(pts[1]);
      nodes.append(pts[2]);
    };
    if (subidx == 4) {
      nodes.append(pts[1]);
      nodes.append(pts[5]);
    };
    if (subidx == 5) {
      nodes.append(pts[2]);
      nodes.append(pts[3]);
    };
    if (subidx == 6) {
      nodes.append(pts[2]);
      nodes.append(pts[6]);
    };
    if (subidx == 7) {
      nodes.append(pts[3]);
      nodes.append(pts[7]);
    };
    if (subidx == 8) {
      nodes.append(pts[4]);
      nodes.append(pts[5]);
    };
    if (subidx == 9) {
      nodes.append(pts[4]);
      nodes.append(pts[7]);
    };
    if (subidx == 10) {
      nodes.append(pts[5]);
      nodes.append(pts[6]);
    };
    if (subidx == 11) {
      nodes.append(pts[6]);
      nodes.append(pts[7]);
    };
  };
};

PolyMesh::PolyMesh(vtkUnstructuredGrid *a_grid, bool dual_mesh)
{
  dbg = false;
  dual = dual_mesh;
  grid = a_grid;
  weight.fill(1.0,grid->GetNumberOfPoints());
  pass1();
  pass2();
  pass3();
};

int PolyMesh::pcIdxNode(vtkIdType id_node)
{
  if (node2pc[id_node] != -1) {
    return node2pc[id_node];
  };
  node2pc[id_node] = id_newpc;
  ++id_newpc;
  return node2pc[id_node];
};

int PolyMesh::pcIdxCell(vtkIdType id_cell)
{
  if (cell2pc[id_cell] != -1) {
    return cell2pc[id_cell];
  };
  cell2pc[id_cell] = id_newpc;
  ++id_newpc;
  return cell2pc[id_cell];
};

void PolyMesh::pass1Tetras()
{
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  face_t F;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_TETRA) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      
      if (dual) {
        
      // edge 0 -> 1
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 0);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 1);
        F.node[3]   = node_t(this, edge, id_cell, 0);
        F.owner     = pcIdxNode(pts[0]);
        F.neighbour = pcIdxNode(pts[1]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // edge 0 -> 2
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 2);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 0);
        F.node[3]   = node_t(this, edge, id_cell, 1);
        F.owner     = pcIdxNode(pts[0]);
        F.neighbour = pcIdxNode(pts[2]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // edge 0 -> 3
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 1);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 2);
        F.node[3]   = node_t(this, edge, id_cell, 2);
        F.owner     = pcIdxNode(pts[0]);
        F.neighbour = pcIdxNode(pts[3]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // edge 1 -> 2
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 0);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 3);
        F.node[3]   = node_t(this, edge, id_cell, 3);
        F.owner     = pcIdxNode(pts[1]);
        F.neighbour = pcIdxNode(pts[2]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // edge 1 -> 3
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 3);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 1);
        F.node[3]   = node_t(this, edge, id_cell, 4);
        F.owner     = pcIdxNode(pts[1]);
        F.neighbour = pcIdxNode(pts[3]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // edge 2 -> 3
        F.node.resize(4);
        F.node[0]   = node_t(this, face, id_cell, 2);
        F.node[1]   = node_t(this, cell, id_cell, 0);
        F.node[2]   = node_t(this, face, id_cell, 3);
        F.node[3]   = node_t(this, edge, id_cell, 5);
        F.owner     = pcIdxNode(pts[2]);
        F.neighbour = pcIdxNode(pts[3]);
        F.bc        = 0;
        F.checkOrientation();
        face_list.append(F);
        
      // boundary face 0
        {
          vtkIdType id_bcell = c2c[id_cell][0];
          if (id_bcell == -1) EG_BUG;
          if (isSurface(id_bcell, grid)) {
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 1);
            F.node[1]   = node_t(this, face, id_cell, 0);
            F.node[2]   = node_t(this, edge, id_cell, 0);
            F.node[3]   = node_t(this, node, pts[0], 0);
            F.owner     = pcIdxNode(pts[0]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 0);
            F.node[1]   = node_t(this, face, id_cell, 0);
            F.node[2]   = node_t(this, edge, id_cell, 3);
            F.node[3]   = node_t(this, node, pts[1], 0);
            F.owner     = pcIdxNode(pts[1]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 3);
            F.node[1]   = node_t(this, face, id_cell, 0);
            F.node[2]   = node_t(this, edge, id_cell, 1);
            F.node[3]   = node_t(this, node, pts[2], 0);
            F.owner     = pcIdxNode(pts[2]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
          };
        };
        
      // boundary face 1
        {
          vtkIdType id_bcell = c2c[id_cell][1];
          if (id_bcell == -1) EG_BUG;
          if (isSurface(id_bcell, grid)) {
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 0);
            F.node[1]   = node_t(this, face, id_cell, 1);
            F.node[2]   = node_t(this, edge, id_cell, 2);
            F.node[3]   = node_t(this, node, pts[0], 0);
            F.owner     = pcIdxNode(pts[0]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 4);
            F.node[1]   = node_t(this, face, id_cell, 1);
            F.node[2]   = node_t(this, edge, id_cell, 0);
            F.node[3]   = node_t(this, node, pts[1], 0);
            F.owner     = pcIdxNode(pts[1]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 2);
            F.node[1]   = node_t(this, face, id_cell, 1);
            F.node[2]   = node_t(this, edge, id_cell, 4);
            F.node[3]   = node_t(this, node, pts[3], 0);
            F.owner     = pcIdxNode(pts[3]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
          };
        };
        
      // boundary face 2
        {
          vtkIdType id_bcell = c2c[id_cell][2];
          if (id_bcell == -1) EG_BUG;
          if (isSurface(id_bcell, grid)) {
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 2);
            F.node[1]   = node_t(this, face, id_cell, 2);
            F.node[2]   = node_t(this, edge, id_cell, 1);
            F.node[3]   = node_t(this, node, pts[0], 0);
            F.owner     = pcIdxNode(pts[0]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 1);
            F.node[1]   = node_t(this, face, id_cell, 2);
            F.node[2]   = node_t(this, edge, id_cell, 5);
            F.node[3]   = node_t(this, node, pts[2], 0);
            F.owner     = pcIdxNode(pts[2]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 5);
            F.node[1]   = node_t(this, face, id_cell, 2);
            F.node[2]   = node_t(this, edge, id_cell, 2);
            F.node[3]   = node_t(this, node, pts[3], 0);
            F.owner     = pcIdxNode(pts[3]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
          };
        };
        
      // boundary face 3
        {
          vtkIdType id_bcell = c2c[id_cell][3];
          if (id_bcell == -1) EG_BUG;
          if (isSurface(id_bcell, grid)) {
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 3);
            F.node[1]   = node_t(this, face, id_cell, 3);
            F.node[2]   = node_t(this, edge, id_cell, 4);
            F.node[3]   = node_t(this, node, pts[1], 0);
            F.owner     = pcIdxNode(pts[1]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 5);
            F.node[1]   = node_t(this, face, id_cell, 3);
            F.node[2]   = node_t(this, edge, id_cell, 3);
            F.node[3]   = node_t(this, node, pts[2], 0);
            F.owner     = pcIdxNode(pts[2]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
            F.node.resize(4);
            F.node[0]   = node_t(this, edge, id_cell, 4);
            F.node[1]   = node_t(this, face, id_cell, 3);
            F.node[2]   = node_t(this, edge, id_cell, 5);
            F.node[3]   = node_t(this, node, pts[3], 0);
            F.owner     = pcIdxNode(pts[3]);
            F.neighbour = -1;
            F.bc        = bc->GetValue(id_bcell);
            face_list.append(F);
            
          };
        };
      } else {
      };
    };
  };
};

void PolyMesh::pass1Prisms()
{
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  face_t F;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_WEDGE) {
      vtkIdType N_pts, *pts, id_ncell;
      bool f0_split = false;
      bool f1_split = false;
      grid->GetCellPoints(id_cell, N_pts, pts);
      
      // face 0
      id_ncell = c2c[id_cell][0];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        F.node.resize(3);
        F.node[0] = node_t(this, node, pts[0], 0);
        F.node[1] = node_t(this, node, pts[1], 0);
        F.node[2] = node_t(this, node, pts[2], 0);
        F.owner   = pcIdxCell(id_cell);
        F.neighbour = -1;
        face_list.append(F);
      } else {
        if (grid->GetCellType(id_ncell) == VTK_WEDGE) {
          F.bc = 0;
          F.node.resize(3);
          F.node[0]   = node_t(this, node, pts[0], 0);
          F.node[1]   = node_t(this, node, pts[1], 0);
          F.node[2]   = node_t(this, node, pts[2], 0);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxCell(id_ncell);
          if (F.owner < F.neighbour) face_list.append(F);
        } else {
          f0_split = true;
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[0], 0);
          F.node[1]   = node_t(this, edge, id_cell, 0);
          F.node[2]   = node_t(this, face, id_cell, 0);
          F.node[3]   = node_t(this, edge, id_cell, 1);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[0]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[1], 0);
          F.node[1]   = node_t(this, edge, id_cell, 3);
          F.node[2]   = node_t(this, face, id_cell, 0);
          F.node[3]   = node_t(this, edge, id_cell, 0);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[1]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[2], 0);
          F.node[1]   = node_t(this, edge, id_cell, 1);
          F.node[2]   = node_t(this, face, id_cell, 0);
          F.node[3]   = node_t(this, edge, id_cell, 3);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[2]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
        };        
      };
      
      // face 1
      id_ncell = c2c[id_cell][1];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        F.node.resize(3);
        F.node[0] = node_t(this, node, pts[3], 0);
        F.node[1] = node_t(this, node, pts[5], 0);
        F.node[2] = node_t(this, node, pts[4], 0);
        F.owner   = pcIdxCell(id_cell);
        F.neighbour = -1;
        face_list.append(F);
      } else {
        if (grid->GetCellType(id_ncell) == VTK_WEDGE) {
          F.bc = 0;
          F.node.resize(3);
          F.node[0]   = node_t(this, node, pts[3], 0);
          F.node[1]   = node_t(this, node, pts[5], 0);
          F.node[2]   = node_t(this, node, pts[4], 0);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxCell(id_ncell);
          if (F.owner < F.neighbour) face_list.append(F);
        } else {
          f1_split = true;
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[3], 0);
          F.node[1]   = node_t(this, edge, id_cell, 7);
          F.node[2]   = node_t(this, face, id_cell, 1);
          F.node[3]   = node_t(this, edge, id_cell, 6);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[3]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[4], 0);
          F.node[1]   = node_t(this, edge, id_cell, 6);
          F.node[2]   = node_t(this, face, id_cell, 1);
          F.node[3]   = node_t(this, edge, id_cell, 8);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[4]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
          
          F.node.resize(4);
          F.node[0]   = node_t(this, node, pts[5], 0);
          F.node[1]   = node_t(this, edge, id_cell, 8);
          F.node[2]   = node_t(this, face, id_cell, 1);
          F.node[3]   = node_t(this, edge, id_cell, 7);
          F.owner     = pcIdxCell(id_cell);
          F.neighbour = pcIdxNode(pts[5]);
          F.bc        = 0;
          F.checkOrientation();
          face_list.append(F);
        };        
      };
      
      // face 2
      id_ncell = c2c[id_cell][2];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.node.resize(4);
      F.owner   = pcIdxCell(id_cell);
      if (!f0_split && !f1_split) {
        F.node.resize(4);
        F.node[0] = node_t(this, node, pts[1], 0);
        F.node[1] = node_t(this, node, pts[0], 0);
        F.node[2] = node_t(this, node, pts[3], 0);
        F.node[3] = node_t(this, node, pts[4], 0);
      };
      if (f0_split && !f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[1], 0);
        F.node[1] = node_t(this, edge, id_cell, 0);
        F.node[2] = node_t(this, node, pts[0], 0);
        F.node[3] = node_t(this, node, pts[3], 0);
        F.node[4] = node_t(this, node, pts[4], 0);
      };
      if (!f0_split && f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[1], 0);
        F.node[1] = node_t(this, node, pts[0], 0);
        F.node[2] = node_t(this, node, pts[3], 0);
        F.node[3] = node_t(this, edge, id_cell, 6);
        F.node[4] = node_t(this, node, pts[4], 0);
      };
      if (f0_split && f1_split) {
        F.node.resize(6);
        F.node[0] = node_t(this, node, pts[1], 0);
        F.node[1] = node_t(this, edge, id_cell, 0);
        F.node[2] = node_t(this, node, pts[0], 0);
        F.node[3] = node_t(this, node, pts[3], 0);
        F.node[4] = node_t(this, edge, id_cell, 6);
        F.node[5] = node_t(this, node, pts[4], 0);
      };
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 3
      id_ncell = c2c[id_cell][3];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner   = pcIdxCell(id_cell);
      if (!f0_split && !f1_split) {
        F.node.resize(4);
        F.node[0] = node_t(this, node, pts[2], 0);
        F.node[1] = node_t(this, node, pts[1], 0);
        F.node[2] = node_t(this, node, pts[4], 0);
        F.node[3] = node_t(this, node, pts[5], 0);
      };
      if (f0_split && !f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[2], 0);
        F.node[1] = node_t(this, edge, id_cell, 3);
        F.node[2] = node_t(this, node, pts[1], 0);
        F.node[3] = node_t(this, node, pts[4], 0);
        F.node[4] = node_t(this, node, pts[5], 0);
      };
      if (!f0_split && f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[2], 0);
        F.node[1] = node_t(this, node, pts[1], 0);
        F.node[2] = node_t(this, node, pts[4], 0);
        F.node[3] = node_t(this, edge, id_cell, 8);
        F.node[4] = node_t(this, node, pts[5], 0);
      };
      if (f0_split && f1_split) {
        F.node.resize(6);
        F.node[0] = node_t(this, node, pts[2], 0);
        F.node[1] = node_t(this, edge, id_cell, 3);
        F.node[2] = node_t(this, node, pts[1], 0);
        F.node[3] = node_t(this, node, pts[4], 0);
        F.node[4] = node_t(this, edge, id_cell, 8);
        F.node[5] = node_t(this, node, pts[5], 0);
      };
      
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 4
      id_ncell = c2c[id_cell][4];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner   = pcIdxCell(id_cell);
      if (!f0_split && !f1_split) {
        F.node.resize(4);
        F.node[0] = node_t(this, node, pts[0], 0);
        F.node[1] = node_t(this, node, pts[2], 0);
        F.node[2] = node_t(this, node, pts[5], 0);
        F.node[3] = node_t(this, node, pts[3], 0);
        F.owner   = pcIdxCell(id_cell);
      };
      if (f0_split && !f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[0], 0);
        F.node[1] = node_t(this, edge, id_cell, 1);
        F.node[2] = node_t(this, node, pts[2], 0);
        F.node[3] = node_t(this, node, pts[5], 0);
        F.node[4] = node_t(this, node, pts[3], 0);
      };
      if (!f0_split && f1_split) {
        F.node.resize(5);
        F.node[0] = node_t(this, node, pts[0], 0);
        F.node[1] = node_t(this, node, pts[2], 0);
        F.node[2] = node_t(this, node, pts[5], 0);
        F.node[3] = node_t(this, edge, id_cell, 7);
        F.node[4] = node_t(this, node, pts[3], 0);
      };
      if (f0_split && f1_split) {
        F.node.resize(6);
        F.node[0] = node_t(this, node, pts[0], 0);
        F.node[1] = node_t(this, edge, id_cell, 1);
        F.node[2] = node_t(this, node, pts[2], 0);
        F.node[3] = node_t(this, node, pts[5], 0);
        F.node[4] = node_t(this, edge, id_cell, 7);
        F.node[5] = node_t(this, node, pts[3], 0);
      };
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
    };
  };
};

void PolyMesh::pass1Hexas()
{
  EG_VTKDCC(vtkIntArray, bc, grid, "cell_code");
  face_t F;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_HEXAHEDRON) {
      vtkIdType N_pts, *pts, id_ncell;
      grid->GetCellPoints(id_cell, N_pts, pts);
      
      // face 0
      id_ncell = c2c[id_cell][0];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[0], 0);
      F.node[1] = node_t(this, node, pts[3], 0);
      F.node[2] = node_t(this, node, pts[2], 0);
      F.node[3] = node_t(this, node, pts[1], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 1
      id_ncell = c2c[id_cell][1];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[4], 0);
      F.node[1] = node_t(this, node, pts[5], 0);
      F.node[2] = node_t(this, node, pts[6], 0);
      F.node[3] = node_t(this, node, pts[7], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 2
      id_ncell = c2c[id_cell][2];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[0], 0);
      F.node[1] = node_t(this, node, pts[1], 0);
      F.node[2] = node_t(this, node, pts[5], 0);
      F.node[3] = node_t(this, node, pts[4], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 3
      id_ncell = c2c[id_cell][3];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[3], 0);
      F.node[1] = node_t(this, node, pts[7], 0);
      F.node[2] = node_t(this, node, pts[6], 0);
      F.node[3] = node_t(this, node, pts[2], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 4
      id_ncell = c2c[id_cell][4];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[0], 0);
      F.node[1] = node_t(this, node, pts[4], 0);
      F.node[2] = node_t(this, node, pts[7], 0);
      F.node[3] = node_t(this, node, pts[3], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
      // face 5
      id_ncell = c2c[id_cell][5];
      F.bc = 0;
      if (id_ncell == -1) EG_BUG;
      if (isSurface(id_ncell, grid)) {
        F.bc = bc->GetValue(id_ncell);
        id_ncell = -1;
      };
      F.owner = pcIdxCell(id_cell);
      F.node.resize(4);
      F.node[0] = node_t(this, node, pts[1], 0);
      F.node[1] = node_t(this, node, pts[2], 0);
      F.node[2] = node_t(this, node, pts[6], 0);
      F.node[3] = node_t(this, node, pts[5], 0);
      if (id_ncell != -1) {
        F.neighbour = pcIdxCell(id_ncell);
        if (F.owner < F.neighbour) face_list.append(F);
      } else {
        F.neighbour = -1;
        face_list.append(F);
      };
      
    };
  };
};

void PolyMesh::sortFaces()
{
  bcs.clear();
  foreach (face_t F, face_list) {
    if (!bcs.contains(F.bc)) {
      bcs.append(F.bc);
    };
  };
  qSort(bcs.begin(),bcs.end());
  QList<face_t> ordered_face_list;
  foreach (int bc, bcs) {
    foreach (face_t F, face_list) {
      if (F.bc == bc) {
        ordered_face_list.append(F);
      };
    };
  };
  faces.resize(ordered_face_list.size());
  qCopy(ordered_face_list.begin(), ordered_face_list.end(), faces.begin());
};

void PolyMesh::pass1()
{   
  cout << "creating polyhedral grid" << endl;
  cell2pc.fill(-1,grid->GetNumberOfCells());
  node2pc.fill(-1,grid->GetNumberOfPoints());
  face_list.clear();
  id_newpc = 0;
  getAllCells(cells, grid);
  createCellToCell(cells, c2c, grid);
  cout << "pass 1 for tetras" << endl;
  pass1Tetras();
  cout << "pass 1 for prisms" << endl;
  pass1Prisms();
  cout << "pass 1 for hexas" << endl;
  pass1Hexas();
  cout << "sorting faces" << endl;
  sortFaces();
};

void PolyMesh::createNodes()
{
  QHash<node_t,int> nodemap;
  {
    int id_node = 0;
    foreach (face_t F, faces) {
      foreach (node_t N, F.node) {
        if (!nodemap.contains(N)) {
          nodemap.insert(N, id_node);
          ++id_node;
        };
      };
    };
  };
  nodes.resize(nodemap.size());
  QVector<bool> nodeset(nodemap.size(), false);
  for (int i_face = 0; i_face < faces.size(); ++i_face) {
    face_t F = faces[i_face];
    for (int i_node = 0; i_node < F.node.size(); ++i_node) {
      node_t N = F.node[i_node];
      int i = nodemap.value(N);
      if (!nodeset[i]) {
        nodeset[i] = true;
        nodes[i] = vec3_t(N.x,N.y,N.z);
      };      
      N.type = node;
      N.idx = i;
      N.subidx = 0;
      F.node[i_node] = N;
    };
    faces[i_face] = F;
  };
};

void PolyMesh::computeNodes()
{
  for (int i_face = 0; i_face < faces.size(); ++i_face) {
    face_t F = faces[i_face];
    for (int i_node = 0; i_node < F.node.size(); ++i_node) {
      node_t N = F.node[i_node];
      if (N.type == node) {
        vec3_t x;
        grid->GetPoint(N.idx, x.data());
        N.x = x[0];
        N.y = x[1];
        N.z = x[2];
      };
      if (N.type == edge) {
        vtkIdType N_pts, *pts;
        vtkIdType type_cell = grid->GetCellType(N.idx);
        grid->GetCellPoints(N.idx, N_pts, pts);
        vec3_t x1, x2;
        double w1=1, w2=1;
        if (type_cell == VTK_TETRA) {
          if (N.subidx == 0) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[1]];
          };
          if (N.subidx == 1) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[2], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[2]];
          };
          if (N.subidx == 2) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[3], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[3]];
          };
          if (N.subidx == 3) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[2], x2.data());
            w1 = weight[pts[1]];
            w2 = weight[pts[2]];
          };
          if (N.subidx == 4) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[3], x2.data());
            w1 = weight[pts[1]];
            w2 = weight[pts[3]];
          };
          if (N.subidx == 5) {
            grid->GetPoint(pts[2], x1.data());
            grid->GetPoint(pts[3], x2.data());
            w1 = weight[pts[2]];
            w2 = weight[pts[3]];
          };
        };
        if (type_cell == VTK_WEDGE) {
          if (N.subidx == 0) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[1]];
          };
          if (N.subidx == 1) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[2], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[2]];
          };
          if (N.subidx == 2) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[3], x2.data());
            w1 = weight[pts[0]];
            w2 = weight[pts[3]];
          };
          if (N.subidx == 3) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[2], x2.data());
            w1 = weight[pts[1]];
            w2 = weight[pts[2]];
          };
          if (N.subidx == 4) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[4], x2.data());
            w1 = weight[pts[1]];
            w2 = weight[pts[4]];
          };
          if (N.subidx == 5) {
            grid->GetPoint(pts[2], x1.data());
            grid->GetPoint(pts[5], x2.data());
            w1 = weight[pts[2]];
            w2 = weight[pts[5]];
          };
          if (N.subidx == 6) {
            grid->GetPoint(pts[3], x1.data());
            grid->GetPoint(pts[4], x2.data());
            w1 = weight[pts[3]];
            w2 = weight[pts[4]];
          };
          if (N.subidx == 7) {
            grid->GetPoint(pts[3], x1.data());
            grid->GetPoint(pts[5], x2.data());
            w1 = weight[pts[3]];
            w2 = weight[pts[5]];
          };
          if (N.subidx == 8) {
            grid->GetPoint(pts[4], x1.data());
            grid->GetPoint(pts[5], x2.data());
            w1 = weight[pts[4]];
            w2 = weight[pts[5]];
          };
        };
        {
          vec3_t x = 1.0/(w1+w2)*(w1*x1+w2*x2);
          N.x = x[0];
          N.y = x[1];
          N.z = x[2];
        };
      };
      if (N.type == face) {
        vtkIdType N_pts, *pts;
        vtkIdType type_cell = grid->GetCellType(N.idx);
        grid->GetCellPoints(N.idx, N_pts, pts);
        if (type_cell == VTK_TETRA) {
          vec3_t x1, x2, x3;
          double w1=1, w2=1, w3=1;
          if (N.subidx == 0) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            grid->GetPoint(pts[2], x3.data());
            w1 = faceW(weight[pts[0]]);
            w2 = faceW(weight[pts[1]]);
            w3 = faceW(weight[pts[2]]);
          };
          if (N.subidx == 1) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            grid->GetPoint(pts[3], x3.data());
            w1 = faceW(weight[pts[0]]);
            w2 = faceW(weight[pts[1]]);
            w3 = faceW(weight[pts[3]]);
          };
          if (N.subidx == 2) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[2], x2.data());
            grid->GetPoint(pts[3], x3.data());
            w1 = faceW(weight[pts[0]]);
            w2 = faceW(weight[pts[2]]);
            w3 = faceW(weight[pts[3]]);
          };
          if (N.subidx == 3) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[2], x2.data());
            grid->GetPoint(pts[3], x3.data());
            w1 = faceW(weight[pts[1]]);
            w2 = faceW(weight[pts[2]]);
            w3 = faceW(weight[pts[3]]);
          };
          {
            vec3_t x = 1.0/(w1+w2+w3)*(w1*x1+w2*x2+w3*x3);
            N.x = x[0];
            N.y = x[1];
            N.z = x[2];
          };
        };
        if (type_cell == VTK_WEDGE) {
          vec3_t x1, x2, x3, x4;
          double w1, w2, w3;
          if (N.subidx == 0) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            grid->GetPoint(pts[2], x3.data());
            w1 = faceW(weight[pts[0]]);
            w2 = faceW(weight[pts[1]]);
            w3 = faceW(weight[pts[2]]);
            {
              vec3_t x = 1.0/(w1+w2+w3)*(w1*x1+w2*x2+w3*x3);
              N.x = x[0];
              N.y = x[1];
              N.z = x[2];
            };
          };
          if (N.subidx == 1) {
            grid->GetPoint(pts[3], x1.data());
            grid->GetPoint(pts[4], x2.data());
            grid->GetPoint(pts[5], x3.data());
            w1 = faceW(weight[pts[3]]);
            w2 = faceW(weight[pts[4]]);
            w3 = faceW(weight[pts[5]]);
            {
              vec3_t x = 1.0/(w1+w2+w3)*(w1*x1+w2*x2+w3*x3);
              N.x = x[0];
              N.y = x[1];
              N.z = x[2];
            };
          };
          if (N.subidx == 2) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[1], x2.data());
            grid->GetPoint(pts[3], x3.data());
            grid->GetPoint(pts[4], x4.data());
            {
              vec3_t x = 0.25*(x1+x2+x3+x4);
              N.x = x[0];
              N.y = x[1];
              N.z = x[2];
            };
          };
          if (N.subidx == 3) {
            grid->GetPoint(pts[1], x1.data());
            grid->GetPoint(pts[2], x2.data());
            grid->GetPoint(pts[4], x3.data());
            grid->GetPoint(pts[5], x4.data());
            {
              vec3_t x = 0.25*(x1+x2+x3+x4);
              N.x = x[0];
              N.y = x[1];
              N.z = x[2];
            };
          };
          if (N.subidx == 4) {
            grid->GetPoint(pts[0], x1.data());
            grid->GetPoint(pts[2], x2.data());
            grid->GetPoint(pts[3], x3.data());
            grid->GetPoint(pts[5], x4.data());
            {
              vec3_t x = 0.25*(x1+x2+x3+x4);
              N.x = x[0];
              N.y = x[1];
              N.z = x[2];
            };
          };
        };
      };
      if (N.type == cell) {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(N.idx, N_pts, pts);
        N.x = 0;
        N.y = 0;
        N.z = 0;
        for (int j = 0; j < N_pts; ++j) {
          vec3_t x;
          grid->GetPoint(pts[j], x.data());
          N.x += x[0];
          N.y += x[1];
          N.z += x[2];
        };
        N.x *= 1.0/N_pts;
        N.y *= 1.0/N_pts;
        N.z *= 1.0/N_pts;
      };
      F.node[i_node] = N;
    };
    faces[i_face] = F;
  };
};

PolyMesh::face_t PolyMesh::combineFaces(QList<face_t> faces)
{
  if (faces.size() < 2) EG_BUG;
  QVector<face_t> src_faces(faces.size());
  QVector<face_t> dst_faces(faces.size());
  qCopy(faces.begin(), faces.end(), src_faces.begin());
  QVector<bool> src_used(src_faces.size(),false);
  int N = src_faces.size()-1;
  face_t F_cur = src_faces[0];
  int i_cur = 0;
  src_used[0] = true;
  //int dir = 0;
  node_t centre, before, after, match_node;
  if (dbg) {
    cout << "\n\ncombining faces:\n";
    for (int i = 0; i < src_faces.size(); ++i) cout << src_faces[i] << endl;
  };
  {
    int i;
    for (i = 0; i < F_cur.node.size(); ++i) {
      if (F_cur[i].type == node) {
        centre = F_cur[i];
        break;
      };
    };
    if (i == F_cur.node.size()) {
      for (i = 0; i < F_cur.node.size(); ++i) {
        if (F_cur[i].type == edge) {
          centre = F_cur[i];
          break;
        };
      };
      if (i == F_cur.node.size()) EG_BUG;
    };
    before = F_cur[i-1];
    after  = F_cur[i+1];
  };
  if (dbg) cout << "check" << endl;
  {
    int i, j = 0;
    match_node = after;
    int loops = 0;
    do {
      for (i = 0; i < src_faces.size(); ++i) {
        if (!src_used[i]) {
          for (j = 0; j < src_faces[i].node.size(); ++j) {
            if (src_faces[i][j] == centre) {
              if (src_faces[i][j-1] == match_node) break;
            };
          };
          if (j != src_faces[i].node.size()) break;
        };
      };
      if (i != src_faces.size()) {
        ++loops;
        F_cur = src_faces[i];
        i_cur = i;
        before = F_cur[j-1];
        after  = F_cur[j+1];
        match_node = after;
        src_used[i] = true;
      }; 
    } while ((i != src_faces.size()) && (loops < src_faces.size()));
  };
  if (dbg) cout << "check" << endl;
  dst_faces[0] = F_cur;
  int i_dst = 1;
  for (int i = 0; i < src_faces.size(); ++i) {
    if (i == i_cur) src_used[i] = true;
    else            src_used[i] = false;
  };
  match_node = before;
  //cout << "centre: " << centre << endl;
  //cout << "match_node: " << match_node << endl;
  //cout << "face 0: " << F_cur << endl;
  if (dbg) cout << "check" << endl;
  while (N > 0) {
    int i, j = 0;
    for (i = 0; i < src_faces.size(); ++i) {
      if (!src_used[i]) {
        for (j = 0; j < src_faces[i].node.size(); ++j) {
          if (src_faces[i][j] == centre) {
            if (src_faces[i][j+1] == match_node) break;
          };
        };
        if (j != src_faces[i].node.size()) break;
      };
    };
    if (i == src_faces.size()) EG_BUG;
    dst_faces[i_dst] = src_faces[i];
    src_used[i] = true;
    //cout << "face " << i_dst << ": " << dst_faces[i_dst] << endl;
    ++i_dst;
    --N;
    F_cur = src_faces[i];
    match_node = F_cur[j-1];
  };
  if (dbg) cout << "check" << endl;
  QList<node_t> new_nodes;
  {
    int i;
    int k = 0;
    for (i = 0; i < dst_faces.size(); ++i) {
      int j;
      for (j = 0; j < dst_faces[i].node.size(); ++j) {
        if (dst_faces[i][j] == centre) break;
      };
      if (j == dst_faces[i].node.size()) EG_BUG;
      //cout << "j=" << j << endl;
      k = j + 1;
      while (!(dst_faces[i][k+1] == dst_faces[i][j])) {
        new_nodes.append(dst_faces[i][k]);
        //cout << "appending " << dst_faces[i][k] << endl;
        k += 1;
      };
    };
    if (!(dst_faces[i-1][k] == new_nodes.first())) {
      new_nodes.append(dst_faces[i-1][k]);
      new_nodes.append(centre);
    };
  };
  if (dbg) cout << "check" << endl;
  face_t F_new;
  F_new.owner = src_faces[0].owner;
  F_new.neighbour = src_faces[0].neighbour;
  F_new.bc = src_faces[0].bc;
  F_new.node.resize(new_nodes.size());
  /*
  if (dir < 0) {
    QList<node_t>::iterator iter = new_nodes.begin();
    int i = new_nodes.size()-1;
    while (iter != new_nodes.end()) {
      F_new.node[i] = *iter;
      ++iter;
      --i;
    };
  } else {
    qCopy(new_nodes.begin(), new_nodes.end(), F_new.node.begin());
  };
  */
  if (dbg) cout << "check" << endl;
  qCopy(new_nodes.begin(), new_nodes.end(), F_new.node.begin());
  //cout << "dir: " << dir << endl;
  //cout << "new face: " << F_new << endl;
  if (dbg) cout << "check" << endl;
  return F_new;
};

void PolyMesh::pass2()
{
  cout << "pass 2" << endl;
  QVector<vec3_t> nodex_save(grid->GetNumberOfPoints());
  QVector<bool> node_changed(grid->GetNumberOfPoints(),false);
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    grid->GetPoint(id_node, nodex_save[id_node].data());
  };
  cout << "computing node coordinates" << endl;
  computeNodes();
  
  cout << "adjusting outer layer prisms" << endl;
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (grid->GetCellType(id_cell) == VTK_WEDGE) {
      bool ok = true;
      if (grid->GetCellType(c2c[id_cell][2]) == VTK_HEXAHEDRON) ok = false;
      if (grid->GetCellType(c2c[id_cell][3]) == VTK_HEXAHEDRON) ok = false;
      if (grid->GetCellType(c2c[id_cell][4]) == VTK_HEXAHEDRON) ok = false;
      if (ok) {
        vtkIdType N_pts, *pts;
        grid->GetCellPoints(id_cell, N_pts, pts);
        vec3_t x[6];
        for (int i = 0; i < 6; ++i) grid->GetPoint(pts[i], x[i].data());
        if (grid->GetCellType(c2c[id_cell][0]) == VTK_TETRA) {
          if (!node_changed[pts[0]]) {
            x[0] = 0.5*(x[0]+x[3]);
            node_changed[pts[0]] = true;
          };
          if (!node_changed[pts[1]]) {
            x[1] = 0.5*(x[1]+x[4]);
            node_changed[pts[1]] = true;
          };
          if (!node_changed[pts[2]]) {
            x[2] = 0.5*(x[2]+x[5]);
            node_changed[pts[2]] = true;
          };
        };
        if (grid->GetCellType(c2c[id_cell][1]) == VTK_TETRA) {
          if (!node_changed[pts[3]]) {
            x[3] = 0.5*(x[0]+x[3]);
            node_changed[pts[3]] = true;
          };
          if (!node_changed[pts[4]]) {
            x[4] = 0.5*(x[1]+x[4]);
            node_changed[pts[4]] = true;
          };
          if (!node_changed[pts[5]]) {
            x[5] = 0.5*(x[2]+x[5]);
            node_changed[pts[5]] = true;
          };
        };
        for (int i = 0; i < 6; ++i) grid->GetPoints()->SetPoint(pts[i], x[i].data());
      };
    };
  };
  for (int i_face = 0; i_face < faces.size(); ++i_face) {
    face_t F = faces[i_face];
    for (int i_node = 0; i_node < F.node.size(); ++i_node) {
      node_t N = F.node[i_node];
      if (N.type == node) {
        vec3_t x;
        grid->GetPoint(N.idx, x.data());
        N.x = x[0];
        N.y = x[1];
        N.z = x[2];
      };
      F.node[i_node] = N;
    };
    faces[i_face] = F;
  };
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    grid->GetPoints()->SetPoint(id_node, nodex_save[id_node].data());
  };
  
  if (faces.size() < 2) return;
  cout << "sorting faces" << endl;
  qStableSort(faces.begin(), faces.end());
  QList<face_t> new_faces;
  int i = 0;
  cout << "combining faces (" << faces.size() << ")" << endl;
  while (i < faces.size()) {
    //cout << i << endl;
    //if (i > 2230000) dbg = true;
    //else dbg = false;
    face_t F = faces[i];
    QList<face_t> combine_faces;
    while ((i < faces.size()) && (faces[i].owner == F.owner) && (faces[i].neighbour == F.neighbour)) {
      combine_faces.append(faces[i]);
      if (dbg) cout << i << endl;
      ++i;
    };
    bool combine = false;
    if (combine_faces.size() > 1) {
      QList<face_t>::iterator j = combine_faces.begin();
      int bc = j->bc;
      combine = true;
      double mincosa = 1.0;
      vec3_t n(0,0,0);
      while (j != combine_faces.end()) {
        if (j->bc != bc) combine = false;
        n += j->normal();
        ++j;
      };
      if (combine) {
        n.normalise();
        j = combine_faces.begin();
        while (j != combine_faces.end()) {
          vec3_t nj = j->normal();
          nj.normalise();
          double cosa = n*nj;
          mincosa = min(cosa,mincosa);
          ++j;
        };
        if ((mincosa < cos(45.0*M_PI/180)) && (bc != 0)) {
          combine = false;
        };
      };
    };
    if (combine) {
      new_faces.append(combineFaces(combine_faces));
    } else {
      foreach (face_t F, combine_faces) {
        new_faces.append(F);
      };
    };
  };
  cout << "sorting faces" << endl;
  faces.resize(new_faces.size());
  qCopy(new_faces.begin(), new_faces.end(), faces.begin());
  qStableSort(faces.begin(), faces.end());
};

void PolyMesh::pass3()
{
  /*
  QVector<vec3_t> n(grid->GetNumberOfPoints(),vec3_t(0,0,0));
  for (int i = 0; i < faces.size(); ++i) {
    for (int j = 0; j < faces[i].node.size(); ++j) {
      node_t N = faces[i][j];
      N.n = vec3_t(0,0,0);
      faces[i].node[j] = N;
    };
  };
  for (int i = 0; i < faces.size(); ++i) {
    vec3_t ni = faces[i].normal();
    for (int j = 0; j < faces[i].node.size(); ++j) {
      node_t N = faces[i][j];
      N.n += ni;
      faces[i].node[j] = N;
    };
  };
  for (int i = 0; i < faces.size(); ++i) {
    for (int j = 0; j < faces[i].node.size(); ++j) {
      node_t N = faces[i][j];
      if (N.type == node) {
        n[N.idx] += N.n;
      };
    };
  };
  for (int i = 0; i < faces.size(); ++i) {
    for (int j = 0; j < faces[i].node.size(); ++j) {
      node_t N = faces[i][j];
      if (N.type == node) {
        N.n = n[N.idx];
      };
      faces[i].node[j] = N;
    };
  };
  for (int i = 0; i < n.size(); ++i) n[i].normalise();
  double w_flat  = 1.0;
  double w_sharp = 5.0;
  double w_max = w_flat;
  for (int i = 0; i < weight.size(); ++i) weight[i] = w_flat;
  for (int i = 0; i < faces.size(); ++i) {
    vec3_t ni = faces[i].normal();
    ni.normalise();
    for (int j = 0; j < faces[i].node.size(); ++j) {
      node_t N = faces[i][j];
      N.n.normalise();
      double f = min(1.0,(1.0-N.n*ni));
      N.w = w_flat + f*(w_sharp - w_flat);
      w_max = max(w_max,N.w);
      if (N.type == node) weight[N.idx] = N.w;
      faces[i].node[j] = N;
    };
  };
  cout << "maximal edge weighting : " << w_max << endl;
  computeNodes();
  */
  cout << "pass 3" << endl;
  createNodes();
  return;
  {
    QVector<bool> is_boundary(nodes.size(),false);
    QVector<bool> del_node(nodes.size(),false);
    QVector<int> face_count(nodes.size(),0);
    foreach (face_t F, faces) {
      foreach (node_t N, F.node) {
        if (N.type != node) EG_BUG;
        if (F.bc != 0) {
          is_boundary[N.idx] = true;
        };
      };
    };
    foreach (face_t F, faces) {
      foreach (node_t N, F.node) {
        if (is_boundary[N.idx]) {
          if (F.bc != 0) {
            ++face_count[N.idx];
          };
        } else {
          ++face_count[N.idx];
        };
      };
    };
    QVector<int> old2new(nodes.size());
    int N = 0;
    {
      int j = 0;
      for (int i = 0; i < nodes.size(); ++i) {
        old2new[i] = j;
        int fc = 2;
        if (face_count[i] <= fc) {
          del_node[i] = true;
          ++N;
        } else {
          ++j;
        };
      };
    };
    cout << "removing " << N << " nodes" << endl;
    for (int i = 0; i < faces.size(); ++i) {
      face_t F = faces[i];
      QList<node_t> new_node;
      foreach (node_t N, F.node) {
        if (!del_node[N.idx]) {
          N.idx = old2new[N.idx];
          new_node.append(N);
        };
      };
      F.node.resize(new_node.size());
      qCopy(new_node.begin(), new_node.end(), F.node.begin());
      faces[i] = F;
    };
    QVector<vec3_t> new_nodes(nodes.size()-N);
    for (int i = 0; i < nodes.size(); ++i) {
      if (!del_node[i]) {
        new_nodes[old2new[i]] = nodes[i];
      };
    };
    nodes.resize(new_nodes.size());
    qCopy(new_nodes.begin(), new_nodes.end(), nodes.begin());
  };
  cout << "done" << endl;
};
