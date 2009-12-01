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
  m_NumIterations          = 5;
  m_NumRelaxations         = 5;
  m_NumBoundaryCorrections = 20;
  m_NumSearch              = 10;
  m_LSearch                = 0.5;
  m_SmoothPrisms           = true;
  m_FOld                   = 0;
  m_FNew                   = 0;
  m_H                      = 1.5;

  getSet("boundary layer", "tetra weighting",                    1.0, m_WTet);
  getSet("boundary layer", "layer height weighting",             1.0, m_WH);
  getSet("boundary layer", "parallel edges weighting",           3.0, m_WPar);
  getSet("boundary layer", "parallel faces weighting",           5.0, m_WN);
  getSet("boundary layer", "similar face area weighting",        5.0, m_WA);
  getSet("boundary layer", "skewness weighting",                 0.0, m_WSkew);
  getSet("boundary layer", "orthogonality weighting",            0.0, m_WOrth);
  getSet("boundary layer", "sharp features on nodes weighting",  8.0, m_WSharp1);
  getSet("boundary layer", "sharp features on nodes exponent",   1.4, m_ESharp1);
  getSet("boundary layer", "sharp features on edges weighting",  3.0, m_WSharp2);
  getSet("boundary layer", "sharp features on edges exponent",   1.3, m_ESharp2);
  getSet("boundary layer", "number of smoothing sub-iterations", 5,   m_NumIterations);
  
  getSet("boundary layer", "angle for sharp features",         45.00, m_CritAngle);

  getSet("boundary layer", "use strict prism checking", true, m_StrictPrismChecking);

  m_CritAngle = GeometryTools::deg2rad(m_CritAngle);

}

void GridSmoother::markNodes()
{
  m_NodeMarked.fill(false,m_Grid->GetNumberOfPoints());
  QVector<bool> new_mark(m_Grid->GetNumberOfPoints());
  for (int i_iterations = 0; i_iterations < 4; ++i_iterations) {
    qCopy(m_NodeMarked.begin(),m_NodeMarked.end(),new_mark.begin());
    for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
      bool mark_cell = false;
      vtkIdType type_cell, N_pts, *pts;
      type_cell = m_Grid->GetCellType(id_cell);
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      if (type_cell == VTK_WEDGE) {
        mark_cell = true;
      } else {
        for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
          if (m_NodeMarked[pts[i_pts]]) {
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
    qCopy(new_mark.begin(),new_mark.end(),m_NodeMarked.begin());
  }
  m_NumMarkedNodes = 0;
  QVector<vtkIdType> nodes = m_Part.getNodes();
  foreach (vtkIdType id_node, nodes) {
    if (id_node < 0) EG_BUG;
    if (id_node > m_Grid->GetNumberOfPoints()) EG_BUG;
    if (m_NodeMarked[id_node]) {
      ++m_NumMarkedNodes;
    }
  }
}

bool GridSmoother::setNewPosition(vtkIdType id_node, vec3_t x_new)
{
  using namespace GeometryTools;

  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  m_Grid->GetPoints()->SetPoint(id_node, x_new.data());
  bool move = true;
  Elements E;

  l2g_t cells = m_Part.getCells();
  l2l_t n2c = m_Part.getN2C();

  foreach (int i_cells, n2c[id_node]) {
    vtkIdType id_cell = cells[i_cells];
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    if (isVolume(id_cell, m_Grid)) {
      if (GeometryTools::cellVA(m_Grid, id_cell) < 0) {
        move = false;
      }
    }

    if (type_cell == VTK_WEDGE && m_StrictPrismChecking) {
      vtkIdType N_pts, *pts;
      vec3_t xtet[4];
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      bool ok = true;
      for (int i = 0; i < 4; ++i) {     // variation
        ok = true;
        for (int j = 0; j < 3; ++j) {   // tetrahedron
          for (int k = 0; k < 4; ++k) { // node
            m_Grid->GetPoint(pts[E.priTet(i,j,k)], xtet[k].data());
          }
          if (GeometryTools::tetraVol(xtet[0], xtet[1], xtet[2], xtet[3]) < 0) {
            ok = false;
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
  if (!move) {
    m_Grid->GetPoints()->SetPoint(id_node, x_old.data());
  }
  return move;
}

void GridSmoother::correctDx(int i_nodes, vec3_t &Dx)
{
  l2g_t nodes = m_Part.getNodes();
  l2l_t n2c = m_Part.getN2C();
  for (int i_boundary_correction = 0; i_boundary_correction < m_NumBoundaryCorrections; ++i_boundary_correction) {
    foreach (vtkIdType id_cell, n2c[i_nodes]) {
      if (isSurface(id_cell, m_Grid)) {
        double A = GeometryTools::cellVA(m_Grid, id_cell);
        if (A > 1e-20) {
          vec3_t n = GeometryTools::cellNormal(m_Grid, id_cell);
          n.normalise();
          Dx -= (n*Dx)*n;
        }
      } else {
        if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
          vtkIdType N_pts, *pts;
          m_Grid->GetCellPoints(id_cell, N_pts, pts);
          vtkIdType id_surf_node = -1;
          if (pts[3] == nodes[i_nodes]) id_surf_node = pts[0];
          if (pts[4] == nodes[i_nodes]) id_surf_node = pts[1];
          if (pts[5] == nodes[i_nodes]) id_surf_node = pts[2];
          if (id_surf_node != -1) {
            vec3_t x0,x1,x2;
            m_Grid->GetPoint(pts[0],x0.data());
            m_Grid->GetPoint(pts[1],x1.data());
            m_Grid->GetPoint(pts[2],x2.data());
            vec3_t a = x1-x0;
            vec3_t b = x2-x0;
            vec3_t c = b-a;
            double L = (a.abs()+b.abs()+c.abs())/3.0;
            vec3_t n = b.cross(a);
            n.normalise();
            vec3_t x_old;
            m_Grid->GetPoint(nodes[i_nodes],x_old.data());
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
  l2g_t nodes = m_Part.getNodes();
  vtkIdType id_node = nodes[i_nodes];
  vec3_t x_old;
  m_Grid->GetPoint(id_node, x_old.data());
  bool moved = false;
  for (int i_relaxation = 0; i_relaxation < m_NumRelaxations; ++i_relaxation) {
    if (setNewPosition(id_node, x_old + Dx)) {
      moved = true;
      break;
    }
    Dx *= 0.5;
  }
  return moved;
}

double GridSmoother::errThickness(double x) 
{
  if (x > 1) x = 2 - x;
  double delta = 0.01;
  double a     = 5.0;
  double err   = 0;
  if      (x < 0) err = 1.0/delta - 1.0/(a+delta);
  else if (x < 1) err = 1/(a*x+delta) - 1.0/(a+delta);
  return err;
}

double GridSmoother::func(vec3_t x)
{
  l2g_t nodes = m_Part.getNodes();
  l2g_t cells = m_Part.getCells();
  l2l_t n2c = m_Part.getN2C();
  l2l_t c2c = m_Part.getC2C();

  vec3_t x_old;
  m_Grid->GetPoint(nodes[m_INodesOpt], x_old.data());
  m_Grid->GetPoints()->SetPoint(nodes[m_INodesOpt], x.data());
  double f       = 0;
  double f13     = 1.0/3.0;
  double f14     = 0.25;
  
  vec3_t n_node(1,0,0);
  QList<vec3_t> n_pri;
  
  foreach (int i_cells, n2c[m_INodesOpt]) {
    vtkIdType id_cell = cells[i_cells];
    if (isVolume(id_cell, m_Grid)) {
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      QVector<vec3_t> xn(N_pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        m_Grid->GetPoint(pts[i_pts],xn[i_pts].data());
      }
      if (type_cell == VTK_TETRA) {
        double L = 0;
        L += (xn[0]-xn[1]).abs();
        L += (xn[0]-xn[2]).abs();
        L += (xn[0]-xn[3]).abs();
        L += (xn[1]-xn[2]).abs();
        L += (xn[1]-xn[3]).abs();
        L += (xn[2]-xn[3]).abs();
        L /= 6;
        double V1 = GeometryTools::cellVA(m_Grid, id_cell, true);
        double V2 = sqrt(1.0/72.0)*L*L*L;
        double e = sqr((V1-V2)/V2);
        f += m_WTet*e;
      }
      if (type_cell == VTK_WEDGE) {
        double L = 0;
        L += (xn[0]-xn[1]).abs();
        L += (xn[0]-xn[2]).abs();
        L += (xn[1]-xn[2]).abs();
        L *= m_H/3.0;
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
        if (nodes[m_INodesOpt] == pts[3]) {
          n_node = xn[3]-xn[0];
          n_pri.append(n_face[1]);
          m_L[nodes[m_INodesOpt]] = L;
        }
        if (nodes[m_INodesOpt] == pts[4]) {
          n_node = xn[4]-xn[1];
          n_pri.append(n_face[1]);
          m_L[nodes[m_INodesOpt]] = L;
        }
        if (nodes[m_INodesOpt] == pts[5]) {
          n_node = xn[5]-xn[2];
          n_pri.append(n_face[1]);
          m_L[nodes[m_INodesOpt]] = L;
        }
        vec3_t v0 = xn[0]-xn[3];
        vec3_t v1 = xn[1]-xn[4];
        vec3_t v2 = xn[2]-xn[5];
        
        double h0 = v0*n_face[0];
        double h1 = v1*n_face[0];
        double h2 = v2*n_face[0];
        if (h0 > 0.5*L) h0 = max(v0.abs(),h0);
        if (h1 > 0.5*L) h1 = max(v1.abs(),h1);
        if (h2 > 0.5*L) h2 = max(v2.abs(),h2);

        {
          
          double e1 = errThickness(h0/L);
          double e2 = errThickness(h1/L);
          double e3 = errThickness(h2/L);
          double e  = max(e1,max(e2,e3));
          //f += w_h*f13*e1;
          //f += w_h*f13*e2;
          //f += w_h*f13*e3;
          //err_pria->SetValue(id_cell,f13*(e1+e2+e3));
          
          f += m_WH*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          v0.normalise();
          v1.normalise();
          v2.normalise();
          double e1 = f13*(1-v0*v1);
          double e2 = f13*(1-v0*v2);
          double e3 = f13*(1-v1*v2);
          f += m_WPar*e1;
          f += m_WPar*e2;
          f += m_WPar*e3;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e = (1+n_face[0]*n_face[1]);
          f += m_WN*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e = sqr((A1-A2)/(A1+A2));
          f += m_WA*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e_skew = 0;
          double e_orth = 0;
          int N = 0;
          vec3_t xc = cellCentre(m_Grid, id_cell);
          for (int i_face = 0; i_face < 5; ++i_face) {
            int i_cells_neigh = c2c[i_cells][i_face];
            if (i_cells_neigh != -1) {
              vtkIdType id_neigh_cell = cells[i_cells_neigh];
              if (isVolume(id_neigh_cell, m_Grid)) {
                vec3_t vc = cellCentre(m_Grid, id_neigh_cell) - xc;
                vec3_t vf = x_face[i_face] - xc;
                vc.normalise();
                vf.normalise();
                e_skew += (1-vc*vf);
                e_orth += (1-vc*n_face[i_face]);
                ++N;
              }
            }
          }
          e_skew /= N;
          e_orth /= N;
          f += m_WSkew*e_skew + m_WOrth*e_orth;
        }
        
        double f_sharp2 = 0;
        for (int j = 2; j <= 4; ++j) {
          vtkIdType id_ncell = c2c[id_cell][j];
          if (id_ncell != -1) {
            if (m_Grid->GetCellType(id_ncell) == VTK_WEDGE) {
              vtkIdType N_pts, *pts;
              m_Grid->GetCellPoints(id_ncell, N_pts, pts);
              QVector<vec3_t> x(3);
              for (int i_pts = 3; i_pts <= 5; ++i_pts) {
                m_Grid->GetPoint(pts[i_pts],x[i_pts-3].data());
              }
              vec3_t n = GeometryTools::triNormal(x[0],x[2],x[1]);
              n.normalise();
              f_sharp2 += pow(fabs(1-n_face[1]*n), m_ESharp2);
            }
          }
        }
        f += m_WSharp2*f_sharp2;
      }
    }
  }
  m_Grid->GetPoints()->SetPoint(nodes[m_INodesOpt], x_old.data());
  n_node.normalise();
  {
    double f_sharp1 = 0;
    foreach (vec3_t n, n_pri) {
      f_sharp1 += pow(fabs(1-n_node*n), m_ESharp1);
    }
    f += m_WSharp1*f_sharp1;
  }
  return f;
}

void GridSmoother::resetStencil()
{
  m_Stencil.clear();
}

void GridSmoother::addToStencil(double C, vec3_t x)
{
  stencil_node_t sn;
  sn.x = x;
  sn.C = C;
  m_Stencil.append(sn);
}

void GridSmoother::operateOptimisation()
{
  markNodes();

  l2g_t nodes = m_Part.getNodes();
  l2g_t cells = m_Part.getCells();
  g2l_t _nodes = m_Part.getLocalNodes();
  l2l_t n2c = m_Part.getN2C();
  l2l_t n2n = m_Part.getN2N();

  EG_VTKDCC(vtkIntArray,    bc,          m_Grid, "cell_code");
  EG_VTKDCN(vtkIntArray,    node_status, m_Grid, "node_status");
  EG_VTKDCN(vtkIntArray,    node_layer,  m_Grid, "node_layer");
  
  m_IdFoot.fill(-1, m_Grid->GetNumberOfPoints());
  m_L.fill(0, m_Grid->GetNumberOfPoints());

  QVector<QSet<int> > n2bc(nodes.size());
  QVector<bool> prism_node(nodes.size(),false);
  foreach (vtkIdType id_cell, cells) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        n2bc[_nodes[pts[i_pts]]].insert(bc->GetValue(id_cell));
      }
    }
    if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
        if (_nodes[pts[i_pts]] != -1) {
          prism_node[_nodes[pts[i_pts]]] = true;
        }
      }
    }
  }
  
  m_FOld = 0;
  m_FMaxOld = 0;
  setPrismWeighting();
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    if (prism_node[i_nodes]) {
      vec3_t x;
      m_INodesOpt = i_nodes;
      m_Grid->GetPoint(nodes[i_nodes], x.data());
      double f = func(x);
      m_FOld += f;
      m_FMaxOld = max(m_FMaxOld, f);
    }
  }
  setAllWeighting();
  
  cout << "\nsmoothing volume mesh (" << m_NumMarkedNodes << " nodes)" << endl;
  for (int i_iterations = 0; i_iterations < m_NumIterations; ++i_iterations) {
    cout << "iteration " << i_iterations+1 << "/" << m_NumIterations << endl;
    int N_blocked  = 0;
    int m_NumSearched = 0;
    int N_illegal  = 0;
    int N1 = 0;
    int N2 = 0;
    int N3 = 0;
    int N4 = 0;
    QTime start = QTime::currentTime();
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      bool smooth = m_NodeMarked[nodes[i_nodes]];
      if (smooth) {
        foreach (int bc, n2bc[i_nodes]) {
          if (m_BoundaryCodes.contains(bc) || (m_BoundaryCodes.size() == 0)) {
            smooth = false;
          }
        }
        foreach (int i_cells, n2c[i_nodes]) {
          vtkIdType type_cell = m_Grid->GetCellType(cells[i_cells]);
          if ((type_cell == VTK_WEDGE) && !m_SmoothPrisms) {
            smooth = false;
          }
        }
      }
      if (smooth) {
        vtkIdType id_node = nodes[i_nodes];
        vtkIdType id_foot = m_IdFoot[id_node];
        bool use_simple = false;
        if (id_foot != -1) {
          if (m_IsSharpNode[id_foot]) {
            use_simple = true;
          }
        }
        if (use_simple) {
          simpleNodeMovement(i_nodes);
          ++N4;
        } else {
          vec3_t xn;
          vec3_t x_old;
          m_Grid->GetPoint(id_node, x_old.data());
          resetStencil();
          bool is_surf = n2bc[i_nodes].size() > 0;
          int N = 0;
          foreach (int j_nodes, n2n[i_nodes]) {
            if (!is_surf || (n2bc[j_nodes].size() > 0)) {
              vtkIdType id_neigh_node = nodes[j_nodes];
              m_Grid->GetPoint(id_neigh_node, xn.data());
              addToStencil(1.0, xn);
              ++N;
            }
          }
          if (N == 0) {
            EG_BUG;
          }
          vec3_t x_new1 = vec3_t(0,0,0);
          m_SumC = 0;
          m_L0 = 0;
          double L_min = 1e99;
          foreach (stencil_node_t sn, m_Stencil) {
            m_SumC += sn.C;
            x_new1 += sn.C*sn.x;
            double L = (sn.x - x_old).abs();
            m_L0 += sn.C*L;
            L_min = min(L_min, L);
          }
          m_L0 /= m_SumC;
          x_new1 *= 1.0/m_SumC;
          vec3_t Dx1 = x_new1 - x_old;
          setDeltas(1e-3*m_L0);
          m_INodesOpt = i_nodes;
          vec3_t Dx2(0,0,0);
          Dx2 = optimise(x_old);
          vec3_t Dx3 = (-10e-4/func(x_old))*grad_f;
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
            vec3_t ex(m_LSearch*m_L0/m_NumSearch, 0, 0);
            vec3_t ey(0,m_LSearch*m_L0/m_NumSearch, 0);
            vec3_t ez(0,0,m_LSearch*m_L0/m_NumSearch);
            vec3_t x_best = x_old;
            double f_min = func(x_old);
            bool found = false;
            bool illegal = false;
            if (!setNewPosition(id_node,x_old)) {
              illegal = true;
            }
            for (int i = -m_NumSearch; i <= m_NumSearch; ++i) {
              for (int j = -m_NumSearch; j <= m_NumSearch; ++j) {
                for (int k = -m_NumSearch; k <= m_NumSearch; ++k) {
                  if ((i != 0) || (j != 0) || (k != 0)) {
                    vec3_t Dx = double(i)*ex + double(j)*ey + double(k)*ez;
                    correctDx(i_nodes, Dx);
                    vec3_t x = x_old + x;
                    if (setNewPosition(id_node, x)) {
                      m_Grid->GetPoints()->SetPoint(id_node, x_old.data());
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
              m_Grid->GetPoints()->SetPoint(id_node, x_best.data());
              ++m_NumSearched;
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
    }
        
    cout << N1 << " type 1 movements (simple)" << endl;
    cout << N2 << " type 2 movements (Newton)" << endl;
    cout << N3 << " type 3 movements (gradient)" << endl;
    //cout << N_blocked << " movements blocked" << endl;
    //cout << m_NumSearched << " movements by search" << endl;
    //cout << N_illegal << " nodes in illegal positions" << endl;
    
    cout << start.secsTo(QTime::currentTime()) << " seconds elapsed" << endl;
    m_FNew = 0;
    m_FMaxNew = 0;
    setPrismWeighting();
    for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
      if (prism_node[i_nodes]) {
        vec3_t x;
        m_INodesOpt = i_nodes;
        m_Grid->GetPoint(nodes[i_nodes], x.data());
        double f = func(x);
        m_FNew += f;
        m_FMaxNew = max(m_FMaxNew, f);
      }
    }
    setAllWeighting();
    cout << "total prism error (old) = " << m_FOld << endl;
    cout << "total prism error (new) = " << m_FNew << endl;
    double f_old     = max(1e-10, m_FOld);
    cout << "total prism improvement = " << 100*(1 - m_FNew/f_old) << "%" << endl;
  }
  cout << "done" << endl;
}

double GridSmoother::improvement()
{
  double f_old = max(1e-10, m_FOld);
  return 1-m_FNew/f_old;
}

void GridSmoother::computeNormals()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  m_NodeNormal.fill(vec3_t(0,0,0), m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    QSet<int> bcs;
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      vtkIdType id_cell = m_Part.n2cGG(id_node, i);
      if (isSurface(id_cell, m_Grid)) {
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
      if (isSurface(id_cell, m_Grid)) {
        int bc = cell_code->GetValue(id_cell);
        if (m_BoundaryCodes.contains(bc)) {
          vtkIdType N_pts, *pts;
          m_Grid->GetCellPoints(id_cell, N_pts, pts);
          vec3_t a, b, c;
          for (int j = 0; j < N_pts; ++j) {
            if (pts[j] == id_node) {
              m_Grid->GetPoint(pts[j], a.data());
              if (j > 0) {
                m_Grid->GetPoint(pts[j-1], b.data());
              } else {
                m_Grid->GetPoint(pts[N_pts-1], b.data());
              }
              if (j < N_pts - 1) {
                m_Grid->GetPoint(pts[j+1], c.data());
              } else {
                m_Grid->GetPoint(pts[0], c.data());
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
  m_IsSharpNode.fill(false, m_Grid->GetNumberOfPoints());
  for (int i = 0; i < m_Part.getNumberOfCells(); ++i) {
    vtkIdType id_cell1 = m_Part.globalCell(i);
    if (isSurface(id_cell1, m_Grid)) {
      if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell1))) {
        vec3_t n1 = cellNormal(m_Grid, id_cell1);
        n1.normalise();
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_cell1, N_pts, pts);
        for (int j = 0; j < m_Part.c2cLSize(i); ++j) {
          vtkIdType id_cell2 = m_Part.c2cLG(i, j);
          if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell2))) {
            vec3_t n2 = cellNormal(m_Grid, id_cell2);
            n2.normalise();
            if (GeometryTools::angle(n1, n2) > m_CritAngle) {
              vtkIdType id_node1 = pts[j];
              vtkIdType id_node2 = pts[0];
              if (j < m_Part.c2cLSize(i) - 1) {
                id_node2 = pts[j+1];
              }
              m_IsSharpNode[id_node1] = true;
              m_IsSharpNode[id_node2] = true;
            }
          }
        }
      }
    }
  }
  m_IsTripleNode.fill(false, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_IsSharpNode[id_node]) {
      QVector<vec3_t> normals;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          if (m_BoundaryCodes.contains(cell_code->GetValue(id_cell))) {
            vec3_t n = cellNormal(m_Grid, id_cell);
            n.normalise();
            normals.push_back(n);
          }
        }
      }
      int N = 0;
      foreach (vec3_t n1, normals) {
        foreach (vec3_t n2, normals) {
          if (GeometryTools::angle(n1, n2) > m_CritAngle) {
            ++N;
          }
        }
      }
    }
  }
}

void GridSmoother::computeFeet()
{
  m_IdFoot.fill(-1, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (m_Grid->GetCellType(id_cell) == VTK_WEDGE) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      m_IdFoot[pts[3]] = pts[0];
      m_IdFoot[pts[4]] = pts[1];
      m_IdFoot[pts[5]] = pts[2];
    }
  }
}

void GridSmoother::simpleNodeMovement(int i_nodes)
{
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired" );
  vtkIdType id_node = m_Part.globalNode(i_nodes);
  vtkIdType id_foot = m_IdFoot[id_node];
  if (id_foot != -1) {
    vec3_t x_foot, x_node;
    m_Grid->GetPoint(id_foot, x_foot.data());
    m_Grid->GetPoint(id_node, x_node.data());
    double L = cl->GetValue(id_foot);
    vec3_t x_new = x_foot + m_RelativeHeight*L*m_NodeNormal[id_foot];
    vec3_t Dx = x_new - x_node;
    correctDx(i_nodes, Dx);
    m_L[id_node] = L;
    moveNode(i_nodes, Dx);
  }
}

void GridSmoother::operateSimple()
{
  cout << "performing simple boundary layer adjustment" << endl;
  computeFeet();
  m_L.fill(0, m_Grid->GetNumberOfPoints());
  l2g_t nodes = m_Part.getNodes();
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    simpleNodeMovement(i_nodes);
  }
}

void GridSmoother::operate()
{
  markNodes();
  computeNormals();
  if (m_SimpleOperation) {
    operateSimple();
  } else {
    operateOptimisation();
  }
}

