// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2011 enGits GmbH                                     +
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
#include "seedsimpleprismaticlayer.h"
#include <vtkIdList.h>

SeedSimplePrismaticLayer::SeedSimplePrismaticLayer()
{
  layer_g  = 0.0;
  layer_dg = 1.0;
}

void SeedSimplePrismaticLayer::setLayerCells(const QVector<vtkIdType> &cells)
{
  layer_cells.resize(cells.size());
  qCopy(cells.begin(), cells.end(), layer_cells.begin());
}

void SeedSimplePrismaticLayer::getLayerCells(QVector<vtkIdType> &cells)
{
  cells.resize(layer_cells.size());
  qCopy(layer_cells.begin(), layer_cells.end(), cells.begin());
}

void SeedSimplePrismaticLayer::prepareLayer()
{
  m_Grid->BuildLinks();
  EG_VTKSP(vtkIdList, nds);
  EG_VTKSP(vtkIdList, cls);
  EG_VTKDCN(vtkIntArray, node_layer,  m_Grid, "node_layer");
  vol_cells.fill(-1, layer_cells.size());
  faces.resize(layer_cells.size());
  N_new_cells  = 0;
  QSet<vtkIdType> new_points;
  for (int i_layer_cell = 0; i_layer_cell < layer_cells.size(); ++i_layer_cell) {
    vtkIdType  id_cell   = layer_cells[i_layer_cell];
    vtkIdType  type_cell = m_Grid->GetCellType(id_cell);
    vtkIdType *pts;
    vtkIdType  Npts;
    m_Grid->GetCellPoints(id_cell, Npts, pts);
    if (type_cell == VTK_TRIANGLE) {
      nds->Reset();
      faces[i_layer_cell].resize(Npts);
      for (int i_pts = 0; i_pts < Npts; ++i_pts) {
        faces[i_layer_cell][Npts - 1 - i_pts] = pts[i_pts];
        nds->InsertNextId(pts[i_pts]);
        new_points.insert(pts[i_pts]);
      }
      m_Grid->GetCellNeighbors(id_cell, nds, cls);
      for (int i_cls = 0; i_cls < cls->GetNumberOfIds(); ++i_cls) {
        if (isVolume(cls->GetId(i_cls), m_Grid)) {
          if (cls->GetId(i_cls) != id_cell) {
            vol_cells[i_layer_cell] = cls->GetId(i_cls);
          }
        }
      }
      N_new_cells  += 1;
    } else if (type_cell == VTK_WEDGE) {
      nds->Reset();
      faces[i_layer_cell].resize(3);
      for (int i_pts = 0; i_pts < 3; ++i_pts) {
        faces[i_layer_cell][2 - i_pts] = pts[i_pts+3];
        //faces[i_layer_cell][i_pts] = pts[i_pts+3];
        nds->InsertNextId(pts[i_pts+3]);
        new_points.insert(pts[i_pts+3]);
      }
      m_Grid->GetCellNeighbors(id_cell, nds, cls);
      for (int i_cls = 0; i_cls < cls->GetNumberOfIds(); ++i_cls) {
        if (isVolume(cls->GetId(i_cls), m_Grid)) {
          if (cls->GetId(i_cls) != id_cell) {
            vol_cells[i_layer_cell] = cls->GetId(i_cls);
          }
        }
      }
      N_new_cells  += 1;
    }
  }
  new_layer = -1;
  old_layer = -1;
  if (faces.size() > 0) {
    old_layer = node_layer->GetValue(faces[0][0]);
    new_layer = old_layer + 1;
  }
  N_new_cells += countBoundaryElements();
  N_new_points = new_points.size();
}

int SeedSimplePrismaticLayer::countBoundaryElements()
{
  int N = 0;
  QVector<QSet<int> > n2f(m_Grid->GetNumberOfPoints());
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      n2f[faces[i_faces][j_faces]].insert(i_faces);
    }
  }
  EG_VTKDCN(vtkIntArray, node_layer,  m_Grid, "node_layer");
  EG_VTKDCN(vtkIntArray, node_status, m_Grid, "node_status");
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      vtkIdType p1 = faces[i_faces][j_faces];
      vtkIdType p2 = faces[i_faces][0];
      if (j_faces < faces[i_faces].size() - 1) {
        p2 = faces[i_faces][j_faces+1];
      }
      bool consider_edge = false;
      if ((node_status->GetValue(p1) & 2) && (node_status->GetValue(p2) & 2)) {
        consider_edge = true;
      }
      if ((node_layer->GetValue(p1) == 0) && (node_layer->GetValue(p2) == 0)) {
        consider_edge = true;
      }
      if (consider_edge) {
        QSet<int> faces_on_edge;
        qcontIntersection(n2f[p1], n2f[p2], faces_on_edge);
        if (faces_on_edge.size() == 0) {
          EG_BUG;
        } else if (faces_on_edge.size() == 1) {
          ++N;
        } else if (faces_on_edge.size() > 2) {
          EG_BUG;
        }
      }
    }
  }
  return N;
}

void SeedSimplePrismaticLayer::createBoundaryElements(vtkUnstructuredGrid *new_grid)
{
  QVector<QSet<int> > n2f(m_Grid->GetNumberOfPoints());
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      n2f[faces[i_faces][j_faces]].insert(i_faces);
    }
  }
  QVector<vtkIdType>  bcells;
  QVector<vtkIdType>  bnodes;
  QVector<int>        _bnodes;
  QVector<QSet<int> > bn2bc;
  getAllSurfaceCells(bcells, new_grid);
  getNodesFromCells(bcells, bnodes, new_grid);
  createNodeMapping(bnodes, _bnodes, new_grid);
  createNodeToCell(bcells, bnodes, _bnodes, bn2bc, new_grid);
  EG_VTKDCC(vtkIntArray, cell_code,   new_grid, "cell_code");
  EG_VTKDCN(vtkIntArray, node_layer,  new_grid, "node_layer");
  EG_VTKDCN(vtkIntArray, node_status, new_grid, "node_status");
  EG_VTKDCC(vtkIntArray, cell_orgdir, new_grid, "cell_orgdir");
  EG_VTKDCC(vtkIntArray, cell_voldir, new_grid, "cell_voldir");
  EG_VTKDCC(vtkIntArray, cell_curdir, new_grid, "cell_curdir");
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      vtkIdType p1 = faces[i_faces][j_faces];
      vtkIdType p2 = faces[i_faces][0];
      if (j_faces < faces[i_faces].size() - 1) {
        p2 = faces[i_faces][j_faces+1];
      }
      bool consider_edge = false;
      if ((node_status->GetValue(p1) & 2) && (node_status->GetValue(p2) & 2)) {
        consider_edge = true;
      }
      if ((node_layer->GetValue(p1) == 0) && (node_layer->GetValue(p2) == 0)) {
        consider_edge = true;
      }
      if (consider_edge) {
        QSet<int> faces_on_edge;
        qcontIntersection(n2f[p1], n2f[p2], faces_on_edge);
        if (faces_on_edge.size() == 0) {
          EG_BUG;
        } else if (faces_on_edge.size() == 1) {
          vtkIdType pts[4] = { p1, p2, old2new[p2], old2new[p1] };
          vtkIdType id_new_bcell = new_grid->InsertNextCell(VTK_QUAD ,4, pts);
          int bc = 9999;
          QSet<int> bc1, bc2, bc3;
          int org_dir = -99;
          int cur_dir = -99;
          int vol_dir = -99;
          if (_bnodes[old2new[p1]] != -1) {
            foreach (int i_bcells, bn2bc[_bnodes[old2new[p1]]]) {
              if (org_dir == -99) {
                org_dir = cell_orgdir->GetValue(bcells[i_bcells]);
              } else if (cell_orgdir->GetValue(bcells[i_bcells]) != org_dir) {
                EG_BUG;
              }
              if (cur_dir == -99) {
                cur_dir = cell_curdir->GetValue(bcells[i_bcells]);
              } else if (cell_curdir->GetValue(bcells[i_bcells]) != cur_dir) {
                EG_BUG;
              }
              if (vol_dir == -99) {
                vol_dir = cell_voldir->GetValue(bcells[i_bcells]);
              } else if (cell_voldir->GetValue(bcells[i_bcells]) != vol_dir) {
                EG_BUG;
              }
              if (!m_BoundaryCodes.contains(cell_code->GetValue(bcells[i_bcells]))) {
                bc1.insert(cell_code->GetValue(bcells[i_bcells]));
              }
            }
          }
          if (_bnodes[old2new[p2]] != -1) {
            foreach (int i_bcells, bn2bc[_bnodes[old2new[p2]]]) {
              if (org_dir == -99) {
                org_dir = cell_orgdir->GetValue(bcells[i_bcells]);
              } else if (cell_orgdir->GetValue(bcells[i_bcells]) != org_dir) {
                EG_BUG;
              }
              if (cur_dir == -99) {
                cur_dir = cell_curdir->GetValue(bcells[i_bcells]);
              } else if (cell_curdir->GetValue(bcells[i_bcells]) != cur_dir) {
                EG_BUG;
              }
              if (vol_dir == -99) {
                vol_dir = cell_voldir->GetValue(bcells[i_bcells]);
              } else if (cell_voldir->GetValue(bcells[i_bcells]) != vol_dir) {
                EG_BUG;
              }
/*              if (!boundary_codes.contains(cell_code->GetValue(bcells[i_bcells]))) {
                bc1.insert(cell_code->GetValue(bcells[i_bcells]));
              }*/
              if (!m_BoundaryCodes.contains(cell_code->GetValue(bcells[i_bcells]))) {
                bc2.insert(cell_code->GetValue(bcells[i_bcells]));
              }
            }
          }
          qcontIntersection(bc1, bc2, bc3);
          if (bc3.size() == 1) {
            bc = *(bc3.begin());
          } else {
            bc = 9999;
          //EG_BUG;
          }
          cell_code->SetValue(id_new_bcell, bc);
          cell_orgdir->SetValue(id_new_bcell, org_dir);
          cell_voldir->SetValue(id_new_bcell, vol_dir);
          cell_curdir->SetValue(id_new_bcell, cur_dir);
          node_status->SetValue(p1, node_status->GetValue(p1) | 2);
          node_status->SetValue(old2new[p1], node_status->GetValue(p1) | 2);
          node_status->SetValue(p2, node_status->GetValue(p2) | 2);
          node_status->SetValue(old2new[p2], node_status->GetValue(p2) | 2);
        } else if (faces_on_edge.size() > 2) {
          EG_BUG;
        }
      }
    }
  }
}

void SeedSimplePrismaticLayer::operate()
{
  cout << "seeding prismatic layer" << endl;
  prepareLayer();
  EG_VTKSP(vtkUnstructuredGrid, new_grid);
  allocateGrid(new_grid, m_Grid->GetNumberOfCells() + N_new_cells, m_Grid->GetNumberOfPoints() + N_new_points);
  vtkIdType id_new_node = 0;
  EG_VTKDCN(vtkIntArray,    node_status_old, m_Grid,   "node_status");
  EG_VTKDCN(vtkIntArray,    node_status_new, new_grid, "node_status");
  EG_VTKDCN(vtkIntArray,    node_layer_old,  m_Grid,   "node_layer");
  EG_VTKDCN(vtkIntArray,    node_layer_new,  new_grid, "node_layer");
  EG_VTKDCC(vtkIntArray,    bc,              m_Grid,   "cell_code");
  EG_VTKDCN(vtkDoubleArray, cl_old,          m_Grid,   "node_meshdensity_desired");
  EG_VTKDCN(vtkDoubleArray, cl_new,          new_grid, "node_meshdensity_desired");

  l2l_t  n2c   = getPartN2C();
  g2l_t _nodes = getPartLocalNodes();
  
  QVector<QSet<int> > n2f(m_Grid->GetNumberOfPoints());
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      n2f[faces[i_faces][j_faces]].insert(i_faces);
    }
  }

  // copy old grid nodes to the new grid
  old2new.resize(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    vec3_t x;
    m_Grid->GetPoints()->GetPoint(id_node, x.data());
    new_grid->GetPoints()->SetPoint(id_new_node, x.data());
    copyNodeData(m_Grid, id_node, new_grid, id_new_node);
    old2new[id_node] = id_new_node;
    ++id_new_node;
  }
  
  QSet<vtkIdType> split_nodes;
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    for (int j_faces = 0; j_faces < faces[i_faces].size(); ++j_faces) {
      split_nodes.insert(faces[i_faces][j_faces]);
    }
  }
  qDebug() << split_nodes.size();
  QVector<bool> is_split_node(m_Grid->GetNumberOfPoints(), false);
  foreach (vtkIdType id_node, split_nodes) {
    vec3_t x;
    m_Grid->GetPoints()->GetPoint(id_node, x.data());
    double h = cl_old->GetValue(id_node);
    
    if (n2f[id_node].size() > 0) {
      vec3_t n(0,0,0);
      double L = 1e99;
      foreach (int i_faces, n2f[id_node]) {
        vec3_t a,b;
        if (faces[i_faces][0] == id_node) {
          m_Grid->GetPoint(faces[i_faces][1],a.data());
          m_Grid->GetPoint(faces[i_faces][2],b.data());
        }
        if (faces[i_faces][1] == id_node) {
          m_Grid->GetPoint(faces[i_faces][2],a.data());
          m_Grid->GetPoint(faces[i_faces][0],b.data());
        }
        if (faces[i_faces][2] == id_node) {
          m_Grid->GetPoint(faces[i_faces][0],a.data());
          m_Grid->GetPoint(faces[i_faces][1],b.data());
        }
        a -= x;
        b -= x;
        L = min(a.abs(),L);
        L = min(b.abs(),L);
        a.normalise();
        b.normalise();
        vec3_t nf = a.cross(b);
        nf.normalise();
        double alpha = acos(a*b);
        n += alpha*nf;
      }
      for (int i_boundary_correction = 0; i_boundary_correction < 20; ++i_boundary_correction) {
        foreach (vtkIdType id_cell, n2c[_nodes[id_node]]) {
          if (isSurface(id_cell, m_Grid)) {
            if (!m_BoundaryCodes.contains(bc->GetValue(id_cell))) {
              double A = GeometryTools::cellVA(m_Grid, id_cell);
              if (A > 1e-60) {
                vec3_t nf = GeometryTools::cellNormal(m_Grid, id_cell);
                nf.normalise();
                n -= (nf*n)*nf;
              }
            }
          }
        }
      }
      n.normalise();
      x += 0.01*L*n;
    } else {
      EG_BUG;
    }
    
    new_grid->GetPoints()->SetPoint(id_new_node, x.data());
    cl_new->SetValue(id_new_node, h);
    old2new[id_node] = id_new_node;
    node_status_new->SetValue(id_new_node, node_status_old->GetValue(id_node));
    node_layer_new->SetValue(id_new_node, node_layer_old->GetValue(id_node) + 1);
    ++id_new_node;
    is_split_node[id_node] = true;
  }
  QVector<bool> needs_correction(m_Grid->GetNumberOfCells(), false);
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if ((type_cell == VTK_TRIANGLE) || (type_cell == VTK_TETRA)) {
      vtkIdType *pts, N_pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      bool split = false;
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        if (is_split_node[pts[i_pts]]) {
          split = true;
          if (node_layer_old->GetValue(pts[i_pts]) != old_layer) {
            EG_BUG;
          }
        }
      }
      if (split) {
        if (type_cell == VTK_TETRA) {
          needs_correction[id_cell] = true;
        } else {
          bool f = false;
          if (old_layer > 0) {
            for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
              if (node_layer_old->GetValue(pts[i_pts]) == old_layer - 1) {
                f = true;
              }
            }
          }
          if (!f) {
            needs_correction[id_cell] = true;
          } else {
            cout << "dodgy face: " << id_cell << "   ";
            for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
              cout << "(" << pts[i_pts] << "," << node_layer_old->GetValue(pts[i_pts]) << ") ";
            }
            cout << endl;
            for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
              if (is_split_node[pts[i_pts]]) {
                vec3_t x;
                m_Grid->GetPoint(pts[i_pts], x.data());
                cout << "split node: " << pts[i_pts] << "  " << x << endl;
              }
            }
          }
        }
      }
    }
  }
  foreach (vtkIdType id_cell, layer_cells) {
    needs_correction[id_cell] = false;
  }
  
  QVector<bool> is_in_vol_cells(m_Grid->GetNumberOfCells(), false);
  for (int i_faces = 0; i_faces < faces.size(); ++i_faces) {
    vtkIdType id_vol_cell = vol_cells[i_faces];
    if (id_vol_cell >= 0) {
      is_in_vol_cells[id_vol_cell] = true;
      vtkIdType type_vol_cell = m_Grid->GetCellType(id_vol_cell);
      if (type_vol_cell == VTK_TETRA) {
        vtkIdType *pts, N_pts, p[6];
        m_Grid->GetCellPoints(id_vol_cell, N_pts, pts);
        p[0] = faces[i_faces][0];
        p[1] = faces[i_faces][1];
        p[2] = faces[i_faces][2];
        p[3] = old2new[p[0]];
        p[4] = old2new[p[1]];
        p[5] = old2new[p[2]];
        {
          vtkIdType pts[6] = { p[2], p[1], p[0], p[5], p[4], p[3] };
          layer_cells[i_faces] = new_grid->InsertNextCell(VTK_WEDGE, 6, pts);
        }
      } else {
        qDebug() << type_vol_cell;
        EG_BUG;
      }
    }
  }
  
//   writeGrid(new_grid, "pre-cellcopy");
  
  // copy old cells to the new grid
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    vtkIdType *pts, N_pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if (needs_correction[id_cell]) {
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        pts[i_pts] = old2new[pts[i_pts]];
      }
    }
    vtkIdType id_new_cell = new_grid->InsertNextCell(type_cell, N_pts, pts);
    copyCellData(m_Grid, id_cell, new_grid, id_new_cell);
  }
  
//   writeGrid(new_grid, "pre-createBoundaryElements");
  
  createBoundaryElements(new_grid);
  UpdateCellIndex(new_grid);
  m_Grid->DeepCopy(new_grid);
  cout << "done." << endl;
}

