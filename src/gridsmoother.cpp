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
#include "gridsmoother.h"
#include "guimainwindow.h"

#include "elements.h"

#include <QTime>

GridSmoother::GridSmoother()
{
  N_iterations = 5;
  N_relaxations = 5;
  N_boundary_corrections = 20;
  N_search = 10;
  L_search = 0.05;
  smooth_prisms = true;
  dbg = false;
  F_old = 0;
  F_new = 0;
  
  getSet("boundary layer", "relative height of boundary layer", 0.75, m_RelativeHeight);
  getSet("boundary layer", "under relaxation",                  1.00, m_UnderRelaxation);
  getSet("boundary layer", "maximal relative edge length",      1.50, m_MaxRelLength);

  getSet("boundary layer", "number of smoothing sub-iterations", 5, N_iterations);

  getSet("boundary layer", "use strict prism checking", true, m_StrictPrismChecking);
  getSet("boundary layer", "write debug file",          true, m_WriteDebugFile);

  getErrSet("boundary layer", "layer height error",      200.0,  2.0, 2.0, 0.10, m_HeightError);
  getErrSet("boundary layer", "edge length error" ,      200.0,  0.1, 2.0, 0.10, m_EdgeLengthError);
  getErrSet("boundary layer", "edge direction error" ,     0.0,  0.0, 2.0, 0.00, m_EdgeDirectionError);
  getErrSet("boundary layer", "tetra error",              50.0,  0.1, 2.0, 0.05, m_TetraError);
  getErrSet("boundary layer", "parallel edges error",      1.0, 10.0, 2.0, 0.00, m_ParallelEdgesError);
  getErrSet("boundary layer", "parallel faces error",      1.0, 20.0, 2.0, 0.00, m_ParallelFacesError);
  getErrSet("boundary layer", "sharp nodes error",         0.0,  0.0, 2.0, 0.00, m_SharpNodesError);
  getErrSet("boundary layer", "sharp edges error",         0.0,  0.0, 2.0, 0.00, m_SharpEdgesError);
  getErrSet("boundary layer", "similar face area error",   0.0,  1.0, 2.0, 0.00, m_SimilarFaceAreaError);
  getErrSet("boundary layer", "feature line error",        0.0,  0.0, 2.0, 0.00, m_FeatureLineError);

  m_CritAngle = GeometryTools::deg2rad(450.0);

}

void GridSmoother::computeNormals()
{
  EG_VTKDCC(vtkIntArray, cell_code, grid, "cell_code");
  m_NodeNormal.fill(vec3_t(0,0,0), grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    QSet<int> bcs;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      if (isSurface(id_cell, grid)) {
        int bc = cell_code->GetValue(id_cell);
        if (m_BoundaryCodes.contains(bc)) {
          bcs.insert(bc);
        }
      }
    }
    int num_bcs = bcs.size();
    QVector<vec3_t> normal(num_bcs, vec3_t(0,0,0));
    QMap<int,int> bcmap;
    int i_bc = 0;
    foreach (int bc, bcs) {
      bcmap[bc] = i_bc;
      ++i_bc;
    }
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      if (isSurface(id_cell, grid)) {
        int bc = cell_code->GetValue(id_cell);
        if (m_BoundaryCodes.contains(bc)) {
          vtkIdType N_pts, *pts;
          grid->GetCellPoints(id_cell, N_pts, pts);
          vec3_t a, b, c;
          for (int j = 0; j < N_pts; ++j) {
            if (pts[j] == id_node) {
              grid->GetPoint(pts[j], a.data());
              if (j > 0) {
                grid->GetPoint(pts[j-1], b.data());
              } else {
                grid->GetPoint(pts[N_pts-1], b.data());
              }
              if (j < N_pts - 1) {
                grid->GetPoint(pts[j+1], c.data());
              } else {
                grid->GetPoint(pts[0], c.data());
              }
            }
          }
          vec3_t u = b - a;
          vec3_t v = c - a;
          double alpha = GeometryTools::angle(u, v);
          vec3_t n = u.cross(v);
          n.normalise();
          normal[bcmap[bc]] += alpha*n;
        }
      }
    }
    for (int i = 0; i < num_bcs; ++i) {
      normal[i].normalise();
    }
    if (num_bcs > 0) {
      if (num_bcs > 1) {
        for (int i = 0; i < num_bcs; ++i) {
          for (int j = i + 1; j < num_bcs; ++j) {
            vec3_t n = normal[i] + normal[j];
            n.normalise();
            m_NodeNormal[id_node] += n;
          }
        }
      } else {
        m_NodeNormal[id_node] = normal[0];
      }
      m_NodeNormal[id_node].normalise();
    }
  }
}

void GridSmoother::markNodes()
{
  node_marked.fill(false,grid->GetNumberOfPoints());
  QVector<bool> new_mark(grid->GetNumberOfPoints());
  for (int i_iterations = 0; i_iterations < 4; ++i_iterations) {
    qCopy(node_marked.begin(),node_marked.end(),new_mark.begin());
    for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
      bool mark_cell = false;
      vtkIdType type_cell, N_pts, *pts;
      type_cell = grid->GetCellType(id_cell);
      grid->GetCellPoints(id_cell, N_pts, pts);
      if (type_cell == VTK_WEDGE) {
        mark_cell = true;
      } else {
        for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
          if (node_marked[pts[i_pts]]) {
            mark_cell = true;
          }
        }
      }
      if (mark_cell) {
        for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
          new_mark[pts[i_pts]] = true;
        }
      }
    }
    qCopy(new_mark.begin(),new_mark.end(),node_marked.begin());
  }
  N_marked_nodes = 0;
  QVector<vtkIdType> nodes = m_Part.getNodes();
  foreach (vtkIdType id_node, nodes) {
    if (id_node < 0) EG_BUG;
    if (id_node > grid->GetNumberOfPoints()) EG_BUG;
    if (node_marked[id_node]) {
      ++N_marked_nodes;
    }
  }
}

bool GridSmoother::setNewPosition(vtkIdType id_node, vec3_t x_new)
{
  using namespace GeometryTools;
  
  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());
  grid->GetPoints()->SetPoint(id_node, x_new.data());
  bool move = true;
  Elements E;
  if (move) {
    l2g_t cells = getPartCells();
    l2l_t n2c   = getPartN2C();
    foreach (int i_cells, n2c[id_node]) {
      vtkIdType id_cell = cells[i_cells];
      vtkIdType type_cell = grid->GetCellType(id_cell);
      if (type_cell == VTK_TETRA) {
        if (GeometryTools::cellVA(grid, id_cell) < 0) {
          move = false;
          //if (dbg) cout << id_node << " : tetra negative" << endl;
        }
      }
      if (type_cell == VTK_WEDGE && m_StrictPrismChecking) {
        vtkIdType N_pts, *pts;
        vec3_t xtet[4];
        grid->GetCellPoints(id_cell, N_pts, pts);
        bool ok = true;
        for (int i = 0; i < 4; ++i) {     // variation
          ok = true;
          for (int j = 0; j < 3; ++j) {   // tetrahedron
            for (int k = 0; k < 4; ++k) { // node
              grid->GetPoint(pts[E.priTet(i,j,k)], xtet[k].data());
            }
            if (GeometryTools::tetraVol(xtet[0], xtet[1], xtet[2], xtet[3]) < 0) {
              ok = false;
              //if (dbg) cout << id_node << " : prism negative" << endl;
            }
          }
          if (ok) {
            break;
          }
        }
        if (!ok) {
          move = false;
        }
      }
    }
  }
  if (!move) {
    grid->GetPoints()->SetPoint(id_node, x_old.data());
  }
  return move;
}

void GridSmoother::correctDx(int i_nodes, vec3_t &Dx)
{
  l2g_t nodes = m_Part.getNodes();
  l2l_t n2c   = m_Part.getN2C();

  for (int i_boundary_correction = 0; i_boundary_correction < N_boundary_corrections; ++i_boundary_correction) {
    foreach (vtkIdType id_cell, n2c[i_nodes]) {
      if (isSurface(id_cell, grid)) {
        double A = GeometryTools::cellVA(grid, id_cell);
        if (A > 1e-20) {
          vec3_t n = GeometryTools::cellNormal(grid, id_cell);
          n.normalise();
          Dx -= (n*Dx)*n;
        }
      } else {
        if (grid->GetCellType(id_cell) == VTK_WEDGE) {
          vtkIdType N_pts, *pts;
          grid->GetCellPoints(id_cell, N_pts, pts);
          vtkIdType id_surf_node = -1;
          if (pts[3] == nodes[i_nodes]) id_surf_node = pts[0];
          if (pts[4] == nodes[i_nodes]) id_surf_node = pts[1];
          if (pts[5] == nodes[i_nodes]) id_surf_node = pts[2];
          if (id_surf_node != -1) {
            vec3_t x0,x1,x2;
            grid->GetPoint(pts[0],x0.data());
            grid->GetPoint(pts[1],x1.data());
            grid->GetPoint(pts[2],x2.data());
            vec3_t a = x1-x0;
            vec3_t b = x2-x0;
            vec3_t c = b-a;
            double L = (a.abs()+b.abs()+c.abs())/3.0;
            vec3_t n = b.cross(a);
            n.normalise();
            vec3_t x_old;
            grid->GetPoint(nodes[i_nodes],x_old.data());
            vec3_t x_new = x_old + Dx - x0;
            if ( (n*x_new) <= 0 ) {
              x_new -= (x_new*n)*n;
              x_new += 1e-4*L*n;
              x_new += x0;
              Dx = x_new - x_old;
            }
          }
        }
      }
    }
  }
}

bool GridSmoother::moveNode(int i_nodes, vec3_t &Dx)
{
  l2g_t nodes = getPartNodes();
  vtkIdType id_node = nodes[i_nodes];
  vec3_t x_old;
  grid->GetPoint(id_node, x_old.data());

  if (m_IdFoot[id_node] != -1) {
    vec3_t x_foot;
    grid->GetPoint(m_IdFoot[id_node], x_foot.data());
    Dx += x_old - x_foot;
    if (Dx.abs() > m_MaxRelLength*m_L[id_node]) {
      Dx.normalise();
      Dx *= m_MaxRelLength*m_L[id_node];
    }
    Dx -= x_old - x_foot;
  }

  bool moved = false;
  for (int i_relaxation = 0; i_relaxation < N_relaxations; ++i_relaxation) {
    if (setNewPosition(id_node, x_old + Dx)) {
      moved = true;
      break;
    }
    Dx *= 0.1;
  }
  return moved;
}

double GridSmoother::func(vec3_t x)
{
  l2g_t nodes = getPartNodes();
  l2l_t n2c   = getPartN2C();
  l2g_t cells = getPartCells();
  l2l_t c2c   = getPartC2C();

  vec3_t x_old;
  grid->GetPoint(nodes[i_nodes_opt], x_old.data());
  grid->GetPoints()->SetPoint(nodes[i_nodes_opt], x.data());
  double f       = 0;
  double f13     = 1.0/3.0;
  double f14     = 0.25;
  
  vec3_t n_node(1,0,0);
  QList<vec3_t> n_pri;

  EG_VTKDCN(vtkDoubleArray, cl, grid, "node_meshdensity_desired" );

  int N_prisms = 0;

  QSet<edge_t> edges;

  foreach (int i_cells, n2c[i_nodes_opt]) {
    vtkIdType id_cell = cells[i_cells];
    if (isVolume(id_cell, grid)) {
      vtkIdType type_cell = grid->GetCellType(id_cell);
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      QVector<vec3_t> xn(N_pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        grid->GetPoint(pts[i_pts],xn[i_pts].data());
      }
      if (type_cell == VTK_TETRA) {
        if (m_TetraError.active()) {
          double L = 0;
          L = max(L, (xn[0]-xn[1]).abs());
          L = max(L, (xn[0]-xn[2]).abs());
          L = max(L, (xn[0]-xn[3]).abs());
          L = max(L, (xn[1]-xn[2]).abs());
          L = max(L, (xn[1]-xn[3]).abs());
          L = max(L, (xn[2]-xn[3]).abs());
          double H = sqrt(2.0/3.0)*L;

          vec3_t n012 = GeometryTools::triNormal(xn[0], xn[1], xn[2]);
          vec3_t n031 = GeometryTools::triNormal(xn[0], xn[3], xn[1]);
          vec3_t n023 = GeometryTools::triNormal(xn[0], xn[2], xn[3]);
          vec3_t n213 = GeometryTools::triNormal(xn[2], xn[1], xn[3]);
          n012.normalise();
          n031.normalise();
          n023.normalise();
          n213.normalise();
          double h = 1e99;
          h = min(h, fabs((xn[0]-xn[3])*n012));
          h = min(h, fabs((xn[0]-xn[2])*n031));
          h = min(h, fabs((xn[0]-xn[1])*n023));
          h = min(h, fabs((xn[0]-xn[3])*n213));
          h /= H;
          f += m_TetraError(h);
        }
      }
      if (type_cell == VTK_WEDGE) {
        ++N_prisms;
        vec3_t a  = xn[2]-xn[0];
        vec3_t b  = xn[1]-xn[0];
        vec3_t c  = xn[5]-xn[3];
        vec3_t d  = xn[4]-xn[3];
        vec3_t n_face[5];
        vec3_t x_face[5];
        n_face[0] = GeometryTools::triNormal(xn[0],xn[1],xn[2]);
        n_face[1] = GeometryTools::triNormal(xn[3],xn[5],xn[4]);
        n_face[2] = GeometryTools::quadNormal(xn[0],xn[3],xn[4],xn[1]);
        n_face[3] = GeometryTools::quadNormal(xn[1],xn[4],xn[5],xn[2]);
        n_face[4] = GeometryTools::quadNormal(xn[0],xn[2],xn[5],xn[3]);
        x_face[0] = f13*(xn[0]+xn[1]+xn[2]);
        x_face[1] = f13*(xn[3]+xn[4]+xn[5]);
        x_face[2] = f14*(xn[0]+xn[3]+xn[4]+xn[1]);
        x_face[3] = f14*(xn[1]+xn[4]+xn[5]+xn[2]);
        x_face[4] = f14*(xn[0]+xn[2]+xn[5]+xn[3]);
        
        double A1 = 0.5*n_face[0].abs();
        double A2 = 0.5*n_face[1].abs();
        for (int i_face = 0; i_face < 5; ++i_face) {
          n_face[i_face].normalise();
        }
        m_IdFoot[pts[3]] = pts[0];
        m_IdFoot[pts[4]] = pts[1];
        m_IdFoot[pts[5]] = pts[2];
        vec3_t v0 = xn[0]-xn[3];
        vec3_t v1 = xn[1]-xn[4];
        vec3_t v2 = xn[2]-xn[5];
        double h0 = v0*n_face[0];
        double h1 = v1*n_face[0];
        double h2 = v2*n_face[0];
        if (m_HeightError.active() || m_EdgeLengthError.active() || m_EdgeDirectionError.active()) {
          double L = 0;
          int i_foot = -1;
          if (nodes[i_nodes_opt] == pts[3]) {
            i_foot = 0;
          }
          if (nodes[i_nodes_opt] == pts[4]) {
            i_foot = 1;
          }
          if (nodes[i_nodes_opt] == pts[5]) {
            i_foot = 2;
          }
          if (i_foot != -1) {
            L = m_RelativeHeight*cl->GetValue(pts[i_foot]);
            n_node = xn[i_foot + 3] - xn[i_foot];
            n_pri.append(n_face[1]);
          }
          double h = 0;
          double l = 0;
          vec3_t ve(0,0,0);
          if (i_foot != -1) {
            m_L[nodes[i_nodes_opt]] = L;
            if (i_foot == 0) {
              h = h0/L;
              l = -(v0*m_NodeNormal[pts[0]])/L;
              ve = v0;
            } else if (i_foot == 1) {
              h = h1/L;
              l = -(v1*m_NodeNormal[pts[1]])/L;
              ve = v1;
            } else if (i_foot == 2) {
              h = h2/L;
              l = -(v2*m_NodeNormal[pts[2]])/L;
              ve = v2;
            }
            ve *= -1;
            ve.normalise();
            f += m_HeightError(h);
            f += m_EdgeLengthError(l);
            f += m_EdgeDirectionError(angleX(m_NodeNormal[pts[i_foot]],ve));
          }
          if (m_ParallelEdgesError.active()) {
            if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
              v0.normalise();
              v1.normalise();
              v2.normalise();
              f += f13*m_ParallelEdgesError(angleX(v0,v1));
              f += f13*m_ParallelEdgesError(angleX(v0,v2));
              f += f13*m_ParallelEdgesError(angleX(v1,v2));
            }
          }
          if (m_ParallelFacesError.active()) {
            if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
              f += m_ParallelFacesError(angleX(-1*n_face[0],n_face[1]));
            }
          }
          if (m_SimilarFaceAreaError.active()) {
            if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
              f += m_SimilarFaceAreaError(fabs(1.0 - 2*(A1-A2)/(A1+A2)));
            }
          }
        }

        if (m_SharpEdgesError.active()) {
          for (int j = 2; j <= 4; ++j) {
            vtkIdType id_ncell = c2c[id_cell][j];
            if (id_ncell != -1) {
              if (grid->GetCellType(id_ncell) == VTK_WEDGE) {
                vtkIdType N_pts, *pts;
                grid->GetCellPoints(id_ncell, N_pts, pts);
                QVector<vec3_t> x(3);
                for (int i_pts = 3; i_pts <= 5; ++i_pts) {
                  grid->GetPoint(pts[i_pts],x[i_pts-3].data());
                }
                vec3_t n = GeometryTools::triNormal(x[0], x[2], x[1]);
                n.normalise();
                f += f13*m_SharpEdgesError(0.5*(1 + n_face[1]*n));
              }
            }
          }
        }
      }
    }
  }
  grid->GetPoints()->SetPoint(nodes[i_nodes_opt], x_old.data());
  n_node.normalise();

  //CHECK THAT THERE ARE NO BOUNDARY NODE NEIGHBOURS!!!!!!!!!!

  if (m_SharpNodesError.active()) {
    if (n_pri.size() > 0) {
      bool apply_error = true;
      for (int i = 0; i < m_Part.n2nLSize(i_nodes_opt); ++i) {
        vtkIdType id_neigh_node = m_Part.n2nLG(i_nodes_opt, i);
        vtkIdType id_foot = m_IdFoot[nodes[i_nodes_opt]];
        int Nbcs = m_Node2BC[id_neigh_node].size();
        if (m_Node2BC[id_neigh_node].size() != 0 && id_neigh_node != m_IdFoot[nodes[i_nodes_opt]]) {
          apply_error = false;
          break;
        }
      }
      if (apply_error) {
        double e = 0;
        foreach (vec3_t n1, n_pri) {
          foreach (vec3_t n2, n_pri) {
            //e = max(e, m_SharpNodesError(angleX(n1, n2)));
            e = max(e, m_SharpNodesError(0.5*(1 + n1*n2)));
          }
        }
        f += e;
      }
    }
  }

  // smooth feature lines
  vtkIdType id_node1 = nodes[i_nodes_opt];
  vtkIdType id_foot1 = m_IdFoot[id_node1];
  if (id_foot1 != -1 && m_FeatureLineError.active()) {
    if (m_Node2BC[id_foot1].size() == 2) {
      QVector<vtkIdType> ids(2, -1);
      int N = 0;
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        vtkIdType id_foot2 = m_IdFoot[id_node2];
        if (id_foot2 != -1) {
          bool found = true;
          foreach (int bc, m_Node2BC[id_foot1]) {
            if (!m_Node2BC[id_foot2].contains(bc)) {
              found = false;
              break;
            }
          }
          if (found) {
            if (N < 2) {
              ids[N] = id_node2;
            }
            ++N;
          }
        }
      }
      if (N == 2) {
        vec3_t a, b;
        grid->GetPoint(ids[0], a.data());
        grid->GetPoint(ids[1], b.data());
        vec3_t u = x - a;
        vec3_t v = b - x;
        f += m_FeatureLineError(0.5*(1 + u*v));
      }
    }
  }

  return f;
}

void GridSmoother::resetStencil()
{
  stencil.clear();
}

void GridSmoother::addToStencil(double C, vec3_t x)
{
  stencil_node_t sn;
  sn.x = x;
  sn.C = C;
  stencil.append(sn);
}

void GridSmoother::operate()
{
  markNodes();
  computeNormals();

  EG_VTKDCC(vtkIntArray, bc,          grid, "cell_code");
  EG_VTKDCN(vtkIntArray, node_status, grid, "node_status");
  EG_VTKDCN(vtkIntArray, node_layer,  grid, "node_layer");

  m_Node2BC.fill(QSet<int>(), grid->GetNumberOfPoints());
  for (vtkIdType id_cell = 0; id_cell < grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, grid)) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        m_Node2BC[pts[i_pts]].insert(bc->GetValue(id_cell));
      }
    }
  }

  m_IdFoot.fill(-1, grid->GetNumberOfPoints());
  m_L.fill(0, grid->GetNumberOfPoints());
  
  
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();
  l2l_t  n2n   = getPartN2N();

  QVector<bool> prism_node(nodes.size(),false);
  foreach (vtkIdType id_cell, cells) {
    if (grid->GetCellType(id_cell) == VTK_WEDGE) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        if (_nodes[pts[i_pts]] != -1) {
          prism_node[_nodes[pts[i_pts]]] = true;
        }
      }
    }
  }
  
  F_old = 0;
  F_max_old = 0;
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    if (prism_node[i_nodes]) {
      vec3_t x;
      i_nodes_opt = i_nodes;
      grid->GetPoint(nodes[i_nodes], x.data());
      double f = func(x);
      F_old += f;
      F_max_old = max(F_max_old,f);
    }
  }
  
  cout << "\nsmoothing volume mesh (" << N_marked_nodes << " nodes)" << endl;
  for (int i_iterations = 0; i_iterations < N_iterations; ++i_iterations) {
    cout << "iteration " << i_iterations+1 << "/" << N_iterations << endl;
    int N_blocked  = 0;
    int N_searched = 0;
    int N_illegal  = 0;
    int N1 = 0;
    int N2 = 0;
    int N3 = 0;
    QTime start = QTime::currentTime();
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      dbg = false;
      bool smooth = node_marked[nodes[i_nodes]];
      if (smooth) {
        foreach (int bc, m_Node2BC[_nodes[i_nodes]]) {
          if (m_BoundaryCodes.contains(bc) || (m_BoundaryCodes.size() ==0)) {
            smooth = false;
          }
        }
        foreach (int i_cells, n2c[i_nodes]) {
          vtkIdType type_cell = grid->GetCellType(cells[i_cells]);
          if ((type_cell == VTK_WEDGE) && !smooth_prisms) {
            smooth = false;
          }
        }
      }
      if (smooth) {
        vtkIdType id_node = nodes[i_nodes];
        vec3_t xn;
        vec3_t x_old;
        grid->GetPoint(id_node, x_old.data());
        resetStencil();
        bool is_surf = m_Node2BC[_nodes[i_nodes]].size() > 0;
        int N = 0;
        foreach (int j_nodes, n2n[i_nodes]) {
          if (!is_surf || (m_Node2BC[_nodes[j_nodes]].size() > 0)) {
            vtkIdType id_neigh_node = nodes[j_nodes];
            grid->GetPoint(id_neigh_node, xn.data());
            addToStencil(1.0, xn);
            ++N;
          }
        }
        if (N == 0) {
          EG_BUG;
        }
        vec3_t x_new1 = vec3_t(0,0,0);
        sum_C = 0;
        L0 = 0;
        double L_min = 1e99;
        foreach (stencil_node_t sn, stencil) {
          sum_C += sn.C;
          x_new1 += sn.C*sn.x;
          double L = (sn.x - x_old).abs();
          L0 += sn.C*L;
          L_min = min(L_min, L);
        }
        L0 /= sum_C;
        x_new1 *= 1.0/sum_C;
        vec3_t Dx1 = x_new1 - x_old;
        setDeltas(1e-6*L0);
        //setDeltas(1e-6);
        i_nodes_opt = i_nodes;
        vec3_t Dx2(0,0,0);
        Dx2 = m_UnderRelaxation*optimise(x_old);
        vec3_t Dx3 = (-1e-4/func(x_old))*grad_f;
        correctDx(i_nodes, Dx1);
        correctDx(i_nodes, Dx2);
        correctDx(i_nodes, Dx3);
        vec3_t Dx = Dx1;
        double f  = func(x_old + Dx1);
        int _N1 = 1;
        int _N2 = 0;
        int _N3 = 0;
        if (f > func(x_old + Dx2)) {
          Dx = Dx2;
          f = func(x_old + Dx2);
          _N1 = 0;
          _N2 = 1;
          _N3 = 0;
        }
        if (f > func(x_old + Dx3)) {
          Dx = Dx3;
          _N1 = 0;
          _N2 = 0;
          _N3 = 1;
        }

        if (!moveNode(i_nodes, Dx)) {
          // search for a better place
          vec3_t x_save = x_old;
          vec3_t ex(L_search*L0/N_search,0,0);
          vec3_t ey(0,L_search*L0/N_search,0);
          vec3_t ez(0,0,L_search*L0/N_search);
          vec3_t x_best = x_old;
          double f_min = func(x_old);
          bool found = false;
          bool illegal = false;
          if (!setNewPosition(id_node,x_old)) {
            illegal = true;
          }
          for (int i = -N_search; i <= N_search; ++i) {
            for (int j = -N_search; j <= N_search; ++j) {
              for (int k = -N_search; k <= N_search; ++k) {
                if ((i != 0) || (j != 0) || (k != 0)) {
                  vec3_t Dx = double(i)*ex + double(j)*ey + double(k)*ez;
                  correctDx(i_nodes, Dx);
                  vec3_t x = x_old + x;
                  if (setNewPosition(id_node, x)) {
                    grid->GetPoints()->SetPoint(id_node, x_old.data());
                    double f = func(x);
                    if (f < f_min) {
                      f_min = f;
                      x_best = x;
                      found = true;
                      illegal = false;
                    }
                  }
                }
              }
            }
          }
          if (found) {
            grid->GetPoints()->SetPoint(id_node, x_best.data());
            ++N_searched;
          } else {
            ++N_blocked;
            if (illegal) {
              ++N_illegal;
            }
          }
        } else {
          N1 += _N1;
          N2 += _N2;
          N3 += _N3;
        }
      }
    }
        
    cout << N1 << " type 1 movements (simple)" << endl;
    cout << N2 << " type 2 movements (Newton)" << endl;
    cout << N3 << " type 3 movements (gradient)" << endl;
    cout << N_searched << " type X movements (search)" << endl;
    cout << N_blocked << " type 0 movements (failure)" << endl;

    cout << start.secsTo(QTime::currentTime()) << " seconds elapsed" << endl;
    F_new = 0;
    F_max_new = 0;

    m_HeightError.reset();
    m_EdgeLengthError.reset();
    m_TetraError.reset();
    m_SharpEdgesError.reset();
    m_SharpNodesError.reset();
    m_ParallelEdgesError.reset();
    m_ParallelFacesError.reset();
    m_SimilarFaceAreaError.reset();
    m_FeatureLineError.reset();
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if (prism_node[i_nodes]) {
        vec3_t x;
        i_nodes_opt = i_nodes;
        grid->GetPoint(nodes[i_nodes], x.data());
        double f = func(x);
        F_new += f;
        F_max_new = max(F_max_new,f);
      }
    }
    cout << "total error (old) = " << F_old << endl;
    cout << "total error (new) = " << F_new << endl;
    double f_old = max(1e-10,F_old);
    cout << "total improvement = " << 100*(1-F_new/f_old) << "%" << endl;
    printMaxErrors();
  }
  if (m_WriteDebugFile) {
    writeDebugFile("gridsmoother");
  }
  cout << "done" << endl;
}

double GridSmoother::improvement()
{
  double f_max_old = max(1e-10,F_max_old);
  double f_old = max(1e-10,F_old);
  double i2 = 1-F_new/f_old;
  return i2;
}

void GridSmoother::printMaxErrors()
{
  cout << "maximal height error = " << m_HeightError.maxError() << endl;
  cout << "maximal edge lenth error = " << m_EdgeLengthError.maxError() << endl;
  cout << "maximal tetra error = " << m_TetraError.maxError() << endl;
  cout << "maximal sharp nodes error = " << m_SharpNodesError.maxError() << endl;
  cout << "maximal sharp edges error = " << m_SharpEdgesError.maxError() << endl;
  cout << "maximal parallel edges error = " << m_ParallelEdgesError.maxError() << endl;
  cout << "maximal parallel faces error = " << m_ParallelFacesError.maxError() << endl;
  cout << "maximal face area error = " << m_SimilarFaceAreaError.maxError() << endl;
  cout << "maximal feature line error = " << m_FeatureLineError.maxError() << endl;
}

void GridSmoother::writeDebugFile(QString file_name)
{
  QVector<vtkIdType> bcells;
  getSurfaceCells(m_BoundaryCodes, bcells, grid);
  MeshPartition bpart(grid);
  bpart.setCells(bcells);
  QVector<vtkIdType> foot2field(bpart.getNumberOfNodes(), -1);
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    vtkIdType id_foot = m_IdFoot[id_node];
    if (id_foot != -1) {
      if (bpart.localNode(id_foot) == -1) {
        EG_BUG;
      }
      foot2field[bpart.localNode(id_foot)] = id_node;
    }
  }
  file_name = GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtk";
  QFile file(file_name);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "# vtk DataFile Version 2.0\n";
  f << "m_NodeNormal\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";
  f << "POINTS " << 2*bpart.getNumberOfNodes() << " float\n";
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vec3_t x;
    vtkIdType id_node = bpart.globalNode(i);
    grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vec3_t x;
    vtkIdType id_node = foot2field[i];
    grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  f << "CELLS " << 2*bpart.getNumberOfCells() << " " << 8*bpart.getNumberOfCells() << "\n";
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType N_pts, *pts;
    if (grid->GetCellType(id_cell) != VTK_TRIANGLE) {
      EG_BUG;
    }
    grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<int> nds1(3), nds2(3);
    for (int j = 0; j < 3; ++j) {
      nds1[j] = bpart.localNode(pts[j]);
      nds2[j] = nds1[j] + bpart.getNumberOfNodes();
    }
    f << "3 " << nds1[0] << " " << nds1[1] << " " << nds1[2] << "\n";
    f << "3 " << nds2[0] << " " << nds2[1] << " " << nds2[2] << "\n";
  }

  f << "CELL_TYPES " << 2*bpart.getNumberOfCells() << "\n";
  for (int i = 0; i < 2*bpart.getNumberOfCells(); ++ i) {
    f << VTK_TRIANGLE << "\n";
  }
  f << "POINT_DATA " << 2*bpart.getNumberOfNodes() << "\n";
  f << "VECTORS N float\n";

  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vtkIdType id_node = bpart.globalNode(i);
    f << m_NodeNormal[id_node][0] << " " << m_NodeNormal[id_node][1] << " " << m_NodeNormal[id_node][2] << "\n";
  }
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    f << "0 0 0\n";
  }
}

