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
  L_search = 0.5;
  smooth_prisms = true;
  dbg = false;
  F_old = 0;
  F_new = 0;
  
  getSet("boundary layer", "tetra weighting", 1.0, w_tet);
  getSet("boundary layer", "layer height weighting", 1.0, w_h);
  getSet("boundary layer", "parallel edges weighting", 3.0, w_par);
  getSet("boundary layer", "parallel faces weighting", 5.0, w_n);
  getSet("boundary layer", "similar face area weighting", 5.0, w_A);
  getSet("boundary layer", "skewness weighting", 0.0, w_skew);
  getSet("boundary layer", "orthogonality weighting", 0.0, w_orth);
  getSet("boundary layer", "sharp features on nodes weighting", 8.0, w_sharp1);
  getSet("boundary layer", "sharp features on nodes exponent", 2.0, e_sharp1);
  getSet("boundary layer", "sharp features on edges weighting", 3.0, w_sharp2);
  getSet("boundary layer", "sharp features on edges exponent", 1.3, e_sharp2);
  getSet("boundary layer", "relative height of boundary layer", 1.5, H);
  getSet("boundary layer", "number of smoothing sub-iterations", 5, N_iterations);
  
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

  QVector<vtkIdType> cells = m_Part.getCells();
  QVector<QVector<int> > n2c = m_Part.getN2C();

  foreach (int i_cells, n2c[id_node]) {
    vtkIdType id_cell = cells[i_cells];
    vtkIdType type_cell = grid->GetCellType(id_cell);
    if (type_cell == VTK_TETRA) {
      if (GeometryTools::cellVA(grid, id_cell) < 0) {
        move = false;
        //if (dbg) cout << id_node << " : tetra negative" << endl;
      }
    }
    if (type_cell == VTK_WEDGE) {
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
  if (!move) {
    grid->GetPoints()->SetPoint(id_node, x_old.data());
  }
  return move;
}

void GridSmoother::correctDx(int i_nodes, vec3_t &Dx)
{
  QVector<vtkIdType> nodes = m_Part.getNodes();
  QVector<QVector<int> > n2c = m_Part.getN2C();

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
  bool moved = false;
  for (int i_relaxation = 0; i_relaxation < N_relaxations; ++i_relaxation) {
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
  /*
  double X[5], Y[5];
  X[0] = 0.0; Y[0] = 10000;
  X[1] = 0.1; Y[1] = 1.0;
  X[2] = 1.0; Y[2] = 0.0;
  X[3] = 2.0; Y[3] = 1.0;
  X[4] = 3.0; Y[4] = 2.0;
  int i = 0;
  while ((i < 3) && (x > X[i+1])) ++i;
  double err = Y[i] + (x-X[i])*(Y[i+1]-Y[i])/(X[i+1]-X[i]);
  if (err < 0) {
    cout << x << ',' << err <<endl;
    cout << i << endl;
    cout << (x-X[i]) << endl;
    cout << (X[i+1]-X[i]) << endl;
    cout << (Y[i+1]-Y[i]) << endl;
    EG_BUG;
  }
  */
  
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
        L += (xn[0]-xn[1]).abs();
        L += (xn[0]-xn[2]).abs();
        L += (xn[0]-xn[3]).abs();
        L += (xn[1]-xn[2]).abs();
        L += (xn[1]-xn[3]).abs();
        L += (xn[2]-xn[3]).abs();
        L /= 6;
        double V1 = GeometryTools::cellVA(grid, id_cell, true);
        double V2 = sqrt(1.0/72.0)*L*L*L;
        double e = sqr((V1-V2)/V2);
        f += w_tet*e;
      }
      if (type_cell == VTK_WEDGE) {
        double L = 0;
        L += (xn[0]-xn[1]).abs();
        L += (xn[0]-xn[2]).abs();
        L += (xn[1]-xn[2]).abs();
        L *= H/3.0;
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
        if (nodes[i_nodes_opt] == pts[3]) {
          n_node = xn[3]-xn[0];
          n_pri.append(n_face[1]);
        }
        if (nodes[i_nodes_opt] == pts[4]) {
          n_node = xn[4]-xn[1];
          n_pri.append(n_face[1]);
        }
        if (nodes[i_nodes_opt] == pts[5]) {
          n_node = xn[5]-xn[2];
          n_pri.append(n_face[1]);
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
          
          f += w_h*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          v0.normalise();
          v1.normalise();
          v2.normalise();
          double e1 = f13*(1-v0*v1);
          double e2 = f13*(1-v0*v2);
          double e3 = f13*(1-v1*v2);
          f += w_par*e1;
          f += w_par*e2;
          f += w_par*e3;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e = (1+n_face[0]*n_face[1]);
          f += w_n*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e = sqr((A1-A2)/(A1+A2));
          f += w_A*e;
        }
        if ((h0 > 0.01*L) && (h1 > 0.01*L) && (h2 > 0.01*L)) {
          double e_skew = 0;
          double e_orth = 0;
          int N = 0;
          vec3_t xc = cellCentre(grid, id_cell);
          for (int i_face = 0; i_face < 5; ++i_face) {
            int i_cells_neigh = c2c[i_cells][i_face];
            if (i_cells_neigh != -1) {
              vtkIdType id_neigh_cell = cells[i_cells_neigh];
              if (isVolume(id_neigh_cell, grid)) {
                vec3_t vc = cellCentre(grid, id_neigh_cell) - xc;
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
          f += w_skew*e_skew + w_orth*e_orth;
        }
        
        double f_sharp2 = 0;
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
              f_sharp2 += pow(fabs(1-n_face[1]*n), e_sharp2);
            }
          }
        }
        f += w_sharp2*f_sharp2;
      }
    }
  }
  grid->GetPoints()->SetPoint(nodes[i_nodes_opt], x_old.data());
  n_node.normalise();
  {
    double f_sharp1 = 0;
    foreach (vec3_t n, n_pri) {
      f_sharp1 += pow(fabs(1-n_node*n), e_sharp1);
    }
    f += w_sharp1*f_sharp1;
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
  setPrismWeighting();
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
  setAllWeighting();
  
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
        foreach (int bc, n2bc[i_nodes]) {
          if (boundary_codes.contains(bc) || (boundary_codes.size() ==0)) {
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
        setDeltas(1e-3*L0);
        i_nodes_opt = i_nodes;
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
    //cout << N_blocked << " movements blocked" << endl;
    //cout << N_searched << " movements by search" << endl;
    //cout << N_illegal << " nodes in illegal positions" << endl;
    
    cout << start.secsTo(QTime::currentTime()) << " seconds elapsed" << endl;
    F_new = 0;
    F_max_new = 0;
    setPrismWeighting();
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
    setAllWeighting();
    cout << "total prism error (old) = " << F_old << endl;
    cout << "total prism error (new) = " << F_new << endl;
    double f_old     = max(1e-10,F_old);
    double f_max_old = max(1e-10,F_max_old);
    cout << "total prism improvement = " << 100*(1-F_new/f_old) << "%" << endl;
    cout << "maximal prism improvement = " << 100*(1-F_max_new/f_max_old) << "%" << endl;
  }
  cout << "done" << endl;
}

double GridSmoother::improvement()
{
  double f_max_old = max(1e-10,F_max_old);
  double i1 = 1-F_max_new/f_max_old;
  double f_old = max(1e-10,F_old);
  double i2 = 1-F_new/f_old;
  return max(i1,i2);
}

