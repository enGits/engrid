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
  
  getSet("boundary layer", "tetra weighting 1",              1000.0,  m_TetraWeighting1);
  getSet("boundary layer", "tetra weighting 2",                 0.1,  m_TetraWeighting2);
  getSet("boundary layer", "tetra switching",                   0.05, m_TetraSwitch);
  getSet("boundary layer", "layer height weighting 1",        100.0,  m_HeightWeighting1);
  getSet("boundary layer", "layer height weighting 2",          2.0,  m_HeightWeighting2);
  getSet("boundary layer", "layer height switching",            0.2,  m_HeightSwitch);
  getSet("boundary layer", "parallel edges weighting",          9.0,  m_ParallelEdgesWeighting);
  getSet("boundary layer", "parallel faces weighting",         15.0,  m_ParallelFacesWeighting);
  getSet("boundary layer", "similar face area weighting",       5.0,  m_SimilarFaceAreaWeighting);
  getSet("boundary layer", "sharp features on nodes weighting", 8.0,  m_SharpNodesWeighting);
  getSet("boundary layer", "sharp features on nodes exponent",  2.0,  m_SharpNodesExponent);
  getSet("boundary layer", "sharp features on edges weighting", 3.0,  m_SharpEdgesWeighting);
  getSet("boundary layer", "sharp features on edges exponent",  1.3,  m_SharpEdgesExponent);
  getSet("boundary layer", "relative height of boundary layer", 0.6,  m_RelativeHeight);
  getSet("boundary layer", "under relaxation",                  1.0,  m_UnderRelaxation);
  getSet("boundary layer", "maximal relative edge length",      1.5,  m_MaxRelLength);

  getSet("boundary layer", "number of smoothing sub-iterations", 5, N_iterations);
  getSet("boundary layer", "use strict prism checking", false, m_StrictPrismChecking);

  m_CritAngle = GeometryTools::deg2rad(450.0);

}

void GridSmoother::computeAngles()
{
  m_MaxAngle.fill(0.0, grid->GetNumberOfPoints());
  m_NodeNormal.fill(vec3_t(0,0,0), grid->GetNumberOfPoints());
  QVector<double> angle_weight(grid->GetNumberOfPoints(), 0.0);
  QVector<bool> angle_set(grid->GetNumberOfPoints(), false);
  for (vtkIdType id_cell1 = 0; id_cell1 < grid->GetNumberOfCells(); ++id_cell1) {
    if (isSurface(id_cell1, grid)) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell1, N_pts, pts);
      if (N_pts == 3) {
        vec3_t n1 = cellNormal(grid, id_cell1);
        for (int i = 0; i < N_pts; ++i) {
          vtkIdType id_node1 = pts[i];
          vtkIdType id_node2 = pts[0];
          vtkIdType id_node3 = pts[1];
          if (i < N_pts - 1) {
            id_node2 = pts[i+1];
            id_node3 = pts[0];
            if (i < N_pts - 2) {
              id_node3 = pts[i+2];
            }
          }
          vtkIdType id_cell2 = m_Part.c2cGG(id_cell1, i);
          if (id_cell2 != -1) {
            vec3_t n2 = cellNormal(grid, id_cell2);
            double alpha = GeometryTools::angle(n1, n2);
            m_MaxAngle[id_node1] = max(alpha, m_MaxAngle[id_node1]);
            m_MaxAngle[id_node2] = max(alpha, m_MaxAngle[id_node2]);
          }
          vec3_t x1, x2, x3;
          grid->GetPoint(id_node1, x1.data());
          grid->GetPoint(id_node2, x2.data());
          grid->GetPoint(id_node3, x3.data());
          vec3_t u = x2 - x1;
          vec3_t v = x3 - x2;
          double alpha = fabs(GeometryTools::angle(x1-x2, x3-x2));
          vec3_t n = v.cross(u);
          n.normalise();
          //if (id_node2 == 18) cout << id_node1 << ',' << id_node2 << ',' << id_node3 << ",   " << n << endl;
          m_NodeNormal[id_node2] += alpha*n;
          angle_weight[id_node2] += alpha;
          angle_set[id_node2] = true;
        }
      }
    }
  }
  for (vtkIdType id_node = 0; id_node < grid->GetNumberOfPoints(); ++id_node) {
    if (angle_set[id_node]) {
      //cout << m_NodeNormal[id_node] << ", " << angle_weight[id_node] << endl;

      m_NodeNormal[id_node] *= 1.0/angle_weight[id_node];
      m_NodeNormal[id_node].normalise();
/*
      if (id_node == 16) {
        cout << m_NodeNormal[id_node] << endl;
        EG_BUG;
      }
*/
    }
  }
}

void GridSmoother::markNodes()
{
  node_marked.fill(false,grid->GetNumberOfPoints());
  QVector<bool> new_mark(grid->GetNumberOfPoints());
  for (int i_iterations = 0; i_iterations < 2; ++i_iterations) {
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

double GridSmoother::errThickness(double x) 
{
  //return fabs(1-x);

  if (x > 1) x = 2 - x;
  const double delta = 0.01;
  const double a     = 5.0;
  const double b     = 1.0;
  const double m0    = -b*a/sqr(a+delta);
  const double err0  = 1.0/delta - 1.0/(a+delta);

  double err   = 0;

  if      (x < 0) err = err0 + x*m0;
  else if (x < 1) err = abs(1/(a*x+delta) - 1.0/(a+delta));
    
  return err/err0;
}

double GridSmoother::errLimit(double x)
{
  const double eps   = 0.1/1.5;
  const double delta = 1.5*eps;
  const double phi   = 0.5/sqr(eps);
  double err = 0;
  if (x < delta) {
    err = phi*sqr(x - delta);
  }
  if (x < .5*eps) {
    err = 1.0 - x/eps;
  }
  return err;
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

  double tetra_error = 0;
  double total_tet_weight = 0;
  bool tets_only = true;
  EG_VTKDCN(vtkDoubleArray, cl, grid, "node_meshdensity_desired" );

  bool base_triangle_found = false;
  bool consider_height_error = false;

  double height1 = 0;
  double height2 = 0;
  //vec3_t n_base(0,0,0);
  vec3_t x_base(0,0,0);
  int N_prisms = 0;
  vtkIdType id_foot = -1;

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
        double e1 = max(0.0, -m_TetraWeighting1*(h - m_TetraSwitch));
        double e2 = fabs(1 - h);
        m_MaxTetError = max(m_MaxTetError, e2);
        tetra_error += max(e1, m_TetraWeighting2*e2);
        total_tet_weight += 1.0;
      }
      if (type_cell == VTK_WEDGE) {
        ++N_prisms;
        tets_only = false;
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
        double L = 0;
        int i_foot = -1;
        if (nodes[i_nodes_opt] == pts[3]) {
          L = m_RelativeHeight*cl->GetValue(pts[0]);
          n_node = xn[3]-xn[0];
          n_pri.append(n_face[1]);
          i_foot = 0;
          id_foot = pts[0];
        }
        if (nodes[i_nodes_opt] == pts[4]) {
          L = m_RelativeHeight*cl->GetValue(pts[1]);
          n_node = xn[4]-xn[1];
          n_pri.append(n_face[1]);
          i_foot = 1;
          id_foot = pts[1];
        }
        if (nodes[i_nodes_opt] == pts[5]) {
          L = m_RelativeHeight*cl->GetValue(pts[2]);
          n_node = xn[5]-xn[2];
          n_pri.append(n_face[1]);
          i_foot = 2;
          id_foot = pts[2];
        }
        m_L[nodes[i_nodes_opt]] = L;
        vec3_t v0 = xn[0]-xn[3];
        vec3_t v1 = xn[1]-xn[4];
        vec3_t v2 = xn[2]-xn[5];
        
        double h0 = v0*n_face[0];
        double h1 = v1*n_face[0];
        double h2 = v2*n_face[0];

        {
          //n_base += n_face[0];
          vec3_t n = n_face[0];
          n.normalise();
          vec3_t v = x - xn[0];
          v -= (v*n)*n;
          vec3_t x0 = xn[0] + v;
          vec3_t x1 = x0 + n;
          vec3_t x2 = x0 - n;
          vec3_t xi, ri;
          if (intersectEdgeAndTriangle(xn[0], xn[1], xn[2], x1, x2, xi, ri, 0.1)) {
            base_triangle_found = true;
          }
        }

        height1 = 1e99;
        if (i_foot != -1) {
          consider_height_error = true;
          if (i_foot == 0) {
            height1 = min(h0/L, height1);
            x_base = xn[0];
            height2 = L;
          } else if (i_foot == 1) {
            height1 = min(h1/L, height1);
            x_base = xn[1];
            height2 = L;
          } else if (i_foot == 2) {
            height1 = min(h2/L, height1);
            x_base = xn[2];
            height2 = L;
          }
        }
        if (m_ParallelEdgesWeighting > 1e-6) {
          if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
            v0.normalise();
            v1.normalise();
            v2.normalise();
            double e1 = f13*(1-v0*v1);
            double e2 = f13*(1-v0*v2);
            double e3 = f13*(1-v1*v2);
            f += m_ParallelEdgesWeighting*e1;
            f += m_ParallelEdgesWeighting*e2;
            f += m_ParallelEdgesWeighting*e3;
            m_MaxParallelEdgesError = max(m_MaxParallelEdgesError, e1 + e2 + e3);
          }
        }
        if (m_ParallelFacesWeighting  > 1e-6) {
          if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
            double e = (1+n_face[0]*n_face[1]);
            f += m_ParallelFacesWeighting*e;
            m_MaxParallelFacesError = max(m_MaxParallelFacesError, e);
          }
        }
        if (m_SimilarFaceAreaWeighting > 1e-6) {
          if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
            double e = sqr((A1-A2)/(A1+A2));
            f += m_SimilarFaceAreaWeighting*e;
            m_MaxFaceAreaError = max(m_MaxFaceAreaError, e);
          }
        }
        if (m_SharpEdgesWeighting > 1e-6) {
          double f_sharp2 = 0;
          int num_sharp2 = 0;
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
                vec3_t n = GeometryTools::triNormal(x[0],x[2],x[1]);
                n.normalise();
                double scal = max(-1.0, min(1.0, n_face[1]*n));
                f_sharp2 += pow(acos(scal)/M_PI, m_SharpEdgesExponent);
                ++num_sharp2;
              }
            }
          }
          if (num_sharp2 > 0) {
            f += m_SharpEdgesWeighting*f_sharp2/num_sharp2;
            m_MaxSharpEdgesError = max(m_MaxSharpEdgesError, f_sharp2/num_sharp2);
          }
        }
      }
    }
  }
  if (!base_triangle_found) {
    //height_amplification = 2;
  }
  double amp = 1.0;
  if (id_foot != -1) {
    if (N_prisms == 0) {
      EG_BUG;
    }
    //n_base.normalise();
    double h = 0;
    //cout << GeometryTools::rad2deg(m_MaxAngle[id_foot]) << endl;
    if (m_NodeNormal[id_foot].abs() < 0.5) EG_BUG;
    double sina = max(1.0/m_MaxRelLength, sin(m_MaxAngle[id_foot]));
    double H = m_RelativeHeight*cl->GetValue(id_foot)/sina;
    h = (x - x_base)*m_NodeNormal[id_foot];
    h /= H;
    double e1 = max(0.0, -m_HeightWeighting1*(h - m_HeightSwitch));
    double e2 = fabs(1 - h);;
    if (fabs(1-h) > m_MaxHeightError) {
      m_MaxHeightError = fabs(1-h);
      m_PosMaxHeightError = x;
    }
    f += amp*max(e1, m_HeightWeighting2*e2);
  }
  grid->GetPoints()->SetPoint(nodes[i_nodes_opt], x_old.data());
  n_node.normalise();
  if (m_SharpNodesWeighting > 1e-6) {
    double f_sharp1 = 0;
    int num_sharp1 = 0;
    foreach (vec3_t n, n_pri) {
      double scal = n_node*n;
      f_sharp1 = max(f_sharp1, pow(acos(scal)/M_PI, m_SharpNodesExponent));
      ++num_sharp1;
    }
    f += m_SharpNodesExponent*f_sharp1;
    m_MaxSharpNodesError = max(m_MaxSharpNodesError, f_sharp1);
  }
  if (tetra_error > 0) {
    tetra_error /= total_tet_weight;
  }
  if (tets_only) {
    f = tetra_error;
  } else {
    f += tetra_error;
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
  computeAngles();
  m_IdFoot.fill(-1, grid->GetNumberOfPoints());
  m_L.fill(0, grid->GetNumberOfPoints());
  
  EG_VTKDCC(vtkIntArray,    bc,          grid, "cell_code");
  EG_VTKDCN(vtkIntArray,    node_status, grid, "node_status");
  EG_VTKDCN(vtkIntArray,    node_layer,  grid, "node_layer");
  
  l2g_t  nodes = getPartNodes();
  g2l_t _nodes = getPartLocalNodes();
  l2g_t  cells = getPartCells();
  l2l_t  n2c   = getPartN2C();
  l2l_t  n2n   = getPartN2N();

  QVector<QSet<int> > n2bc(nodes.size());
  QVector<bool> prism_node(nodes.size(),false);
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, grid)) {
      vtkIdType N_pts, *pts;
      grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        n2bc[_nodes[pts[i_pts]]].insert(bc->GetValue(id_cell));
      }
    }
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
  m_MaxHeightError = 0;
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
  m_MaxHeightError = 0;
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
        foreach (int bc, n2bc[i_nodes]) {
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
        bool is_surf = n2bc[i_nodes].size() > 0;
        int N = 0;
        foreach (int j_nodes, n2n[i_nodes]) {
          if (!is_surf || (n2bc[j_nodes].size() > 0)) {
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
    //setPrismWeighting();
    m_MaxHeightError = 0;
    m_MaxTetError = 0;
    m_MaxSharpNodesError = 0;
    m_MaxSharpEdgesError = 0;
    m_MaxParallelEdgesError = 0;
    m_MaxParallelFacesError = 0;
    m_MaxFaceAreaError = 0;
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
    //setAllWeighting();
    cout << "total prism error (old) = " << F_old << endl;
    cout << "total prism error (new) = " << F_new << endl;
    double f_old = max(1e-10,F_old);
    cout << "total prism improvement = " << 100*(1-F_new/f_old) << "%" << endl;
    cout << "maximal height error = " << m_MaxHeightError << endl;
    cout << "maximal tetra error = " << m_MaxTetError << endl;
    cout << "maximal sharp nodes error = " << m_MaxSharpNodesError << endl;
    cout << "maximal sharp edges error = " << m_MaxSharpEdgesError << endl;
    cout << "maximal parallel edges error = " << m_MaxParallelEdgesError << endl;
    cout << "maximal parallel faces error = " << m_MaxParallelFacesError << endl;
    cout << "maximal face area error = " << m_MaxFaceAreaError << endl;
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

