// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of enGrid.                                         +
// +                                                                      +
// + Copyright 2008-2012 enGits GmbH                                     +
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
#include "optimisenormalvector.h"
#include "pointfinder.h"

#include <QTime>

GridSmoother::GridSmoother()
{
  m_NumIterations          = 5;
  m_NumRelaxations         = 5;
  m_NumBoundaryCorrections = 50;
  m_DesiredStretching      = 1.2;
  m_FirstCall = true;

  getSet("boundary layer", "number of smoothing sub-iterations",       5,     m_NumIterations);
  getSet("boundary layer", "use strict prism checking",                false, m_StrictPrismChecking);
  getSet("boundary layer", "number of normal vector relax iterations", 10,    m_NumNormalRelaxations);
  getSet("boundary layer", "number of layer height relax iterations",  3,     m_NumHeightRelaxations);
  getSet("boundary layer", "radar angle",                              45,    m_RadarAngle);
  getSet("boundary layer", "maximal layer height in gaps",             0.2,   m_MaxHeightInGaps);

  //m_CritAngle = GeometryTools::deg2rad(m_CritAngle);
}

void GridSmoother::markNodes()
{
  m_NodeMarked.fill(false,m_Grid->GetNumberOfPoints());
  QVector<bool> new_mark(m_Grid->GetNumberOfPoints());
  for (int i_iterations = 0; i_iterations < 1; ++i_iterations) {
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
  QSet<int> free_bcs = m_BoundaryCodes + m_LayerAdjacentBoundaryCodes;
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      if (!free_bcs.contains(cell_code->GetValue(id_cell))) {
        vtkIdType N_pts, *pts;
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        for (int i_pts = 0; i_pts < N_pts; ++i_pts) {
          m_NodeMarked[pts[i_pts]] = false;
        }
      }
    }
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

  if (move) {
    Elements E;

    l2g_t cells = m_Part.getCells();
    l2l_t n2c = m_Part.getN2C();

    foreach (int i_cells, n2c[id_node]) {
      vtkIdType id_cell = cells[i_cells];
      vtkIdType type_cell = m_Grid->GetCellType(id_cell);

      if (type_cell == VTK_TRIANGLE) {
        vtkIdType N_pts, *pts;
        vec3_t x[3];
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        for (int i = 0; i < 3; ++i) {
          m_Grid->GetPoint(pts[i], x[i].data());
        }
        double L_max = 0;
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            if (i != j) {
              L_max = max(L_max, (x[i]-x[j]).abs());
            }
          }
        }
        double A = GeometryTools::cellVA(m_Grid, id_cell);
        if (A < 1e-3*L_max*L_max) {
          move = false;
        } else if (m_NodeNormal[id_node].abs() > 0.1) {
          vec3_t n = GeometryTools::triNormal(x[0], x[1], x[2]);
          n.normalise();
          if (n*m_NodeNormal[id_node] < 0.1) {
            move = false;
          }
        }
      }

      if (type_cell == VTK_TETRA) {
        vtkIdType N_pts, *pts;
        vec3_t x[4];
        m_Grid->GetCellPoints(id_cell, N_pts, pts);
        for (int i = 0; i < 4; ++i) {
          m_Grid->GetPoint(pts[i], x[i].data());
        }
        double L_max = 0;
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            if (i != j) {
              L_max = max(L_max, (x[i]-x[j]).abs());
            }
          }
        }
        if (GeometryTools::cellVA(m_Grid, id_cell) < 1e-3*L_max*L_max*L_max) {
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
            double V = GeometryTools::tetraVol(xtet[0], xtet[1], xtet[2], xtet[3]);
            if (V <= 0) {
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
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  vec3_t x_old;
  m_Grid->GetPoint(nodes[i_nodes],x_old.data());
  for (int i_boundary_correction = 0; i_boundary_correction < m_NumBoundaryCorrections; ++i_boundary_correction) {
    foreach (vtkIdType id_cell, n2c[i_nodes]) {
      if (isSurface(id_cell, m_Grid)) {
        int bc = cell_code->GetValue(id_cell);
        vec3_t x_new = x_old + Dx;
        x_new = GuiMainWindow::pointer()->getSurfProj(bc)->projectRestricted(x_new, nodes[i_nodes]);
        Dx = x_new - x_old;
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
  m_CollisionDetected = false;
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

void GridSmoother::computeNormals()
{
  EG_VTKDCC(vtkIntArray, cell_code, m_Grid, "cell_code");
  m_NodeNormal.fill(vec3_t(0,0,0), m_Grid->GetNumberOfPoints());
  QVector<int> num_bcs(m_Grid->GetNumberOfPoints());
  QVector<OptimiseNormalVector> n_opt(m_Grid->GetNumberOfPoints());
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
    num_bcs[id_node] = bcs.size();
    QVector<vec3_t> normal(num_bcs[id_node], vec3_t(0,0,0));
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
        if (m_BoundaryCodes.contains(bc)) {
          normal[bcmap[bc]] += alpha*n;
          n_opt[id_node].addFace(n);
        } else {
          n_opt[id_node].addConstraint(n);
        }
      }
    }
    for (int i = 0; i < num_bcs[id_node]; ++i) {
      normal[i].normalise();
    }
    if (num_bcs[id_node] > 0) {
      if (num_bcs[id_node] > 1) {
        if (num_bcs[id_node] == 3) {
          for (int i = 0; i < num_bcs[id_node]; ++i) {
            for (int j = i + 1; j < num_bcs[id_node]; ++j) {
              vec3_t n = normal[i] + normal[j];
              n.normalise();
              m_NodeNormal[id_node] += n;
            }
          }
        } else {
          for (int i = 0; i < num_bcs[id_node]; ++i) {
            m_NodeNormal[id_node] += normal[i];
          }
        }
      } else {
        m_NodeNormal[id_node] = normal[0];
      }
      m_NodeNormal[id_node].normalise();
      m_NodeNormal[id_node] = n_opt[id_node](m_NodeNormal[id_node]);
      m_NodeNormal[id_node].normalise();
    }
  }

  m_GeoNormal = m_NodeNormal;
  relaxNormalVectors();

  /*
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_NodeNormal[id_node].abs() < 0.1) {
      vec3_t n(0,0,0);
      int N = 0;
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          n += GeometryTools::cellNormal(m_Grid, id_cell);
          ++N;
        }
      }
      if (N) {
        n.normalise();
        m_NodeNormal[id_node] = n;
      }
    }
    if (num_bcs[id_node] > 1) {
      m_NodeNormal[id_node] = n_opt[id_node](m_NodeNormal[id_node]);
    }
  }
  */

}

void GridSmoother::relaxNormalVectors()
{
  m_SurfNode.fill(false, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
      if (isSurface(m_Part.n2cGG(id_node,i), m_Grid)) {
        m_SurfNode[id_node] = true;
        break;
      }
    }
  }
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  QVector<QSet<int> > n2bc(m_Grid->GetNumberOfPoints());
  for (vtkIdType id_cell = 0; id_cell < m_Grid->GetNumberOfCells(); ++id_cell) {
    if (isSurface(id_cell, m_Grid)) {
      vtkIdType N_pts, *pts;
      m_Grid->GetCellPoints(id_cell, N_pts, pts);
      for (int i = 0; i < N_pts; ++i) {
        if (m_SurfNode[pts[i]] && m_BoundaryCodes.contains(bc->GetValue(id_cell))) {
          n2bc[pts[i]].insert(bc->GetValue(id_cell));
        }
      }
    }
  }
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (n2bc[id_node].size() == 0) {
      m_SurfNode[id_node] = false;
    }
  }
  for (int iter = 0; iter < m_NumNormalRelaxations; ++iter) {
    QVector<vec3_t> n_new(m_NodeNormal.size(), vec3_t(0,0,0));
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_SurfNode[id_node] ) {
        QList<vtkIdType> snap_points;
        for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node,i);          
          bool append = true;
          foreach (int bc, n2bc[id_node]) {
            if (!n2bc[id_neigh].contains(bc)) {
              append = false;
              break;
            }
          }
          if (append) {
            if (!m_SurfNode[id_neigh]) {
              EG_BUG;
            }
            snap_points.append(id_neigh);
          }
        }
        if (snap_points.size() > 0) {
          n_new[id_node] = vec3_t(0,0,0);
          foreach (vtkIdType id_snap, snap_points) {
            n_new[id_node] += m_NodeNormal[id_snap];
          }
          n_new[id_node].normalise();
        } else {
          n_new[id_node] = m_NodeNormal[id_node];
        }
        if (n_new[id_node].abs() < 0.1) {
          EG_BUG;
        }
      }
    }
    m_NodeNormal = n_new;
    correctNormalVectors();
  }
}

void GridSmoother::getRules()
{
  QString rules_text = GuiMainWindow::pointer()->getXmlSection("engrid/blayer/rules");
  QStringList rules = rules_text.split(";", QString::SkipEmptyParts);
  foreach (QString rule, rules) {
    rule = rule.trimmed();
    QStringList parts = rule.split("=");
    if (parts.count() > 1) {
      rule_t rule;
      QString left = parts[0].trimmed();
      rule.h = parts[1].trimmed().toDouble();
      QStringList rows = left.split("<OR>");
      foreach (QString row, rows) {
        QStringList cols = row.split("<AND>");
        foreach (QString col, cols) {
          rule.bcs.insert(col.toInt());
        }
      }
    }
  }
}

void GridSmoother::computeDesiredHeights()
{
  // first pass (intial height)
  EG_VTKDCN(vtkDoubleArray, cl, m_Grid, "node_meshdensity_desired" );
  m_Height.fill(0, m_Grid->GetNumberOfPoints());
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_SurfNode[id_node]) {
      m_Height[id_node] = cl->GetValue(id_node);
      // if undefined: compute height from surrounding edges
      if (m_Height[id_node] < 1e-99) {
        m_Height[id_node] = 0;
        int N = 0;
        vec3_t x;
        m_Grid->GetPoint(id_node, x.data());
        for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          vtkIdType id_neigh_node = m_Part.n2nGG(id_node, i);
          if (m_SurfNode[id_neigh_node]) {
            ++N;
            vec3_t xn;
            m_Grid->GetPoint(id_neigh_node, xn.data());
            m_Height[id_node] += (x - xn).abs();
          }
        }
        if (N == 0) {
          EG_BUG;
        }
        m_Height[id_node] /= N;
      }
    }
  }
  QVector<double> Dh_max = m_Height;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    Dh_max[id_node] *= m_FarRatio;
  }

  getRules();
  // second pass (correct with absolute height if required)
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_SurfNode[id_node]) {
      m_Height[id_node] = m_Blending*m_AbsoluteHeight + (1.0-m_Blending)*m_RelativeHeight*m_Height[id_node];
      double h_rule = 1e99;
      foreach (rule_t rule, m_Rules) {
        bool apply_rule = true;
        foreach (int bc, rule.bcs) {
          if (!m_Part.hasBC(id_node, bc)) {
            apply_rule = false;
            break;
          }
        }
        if (apply_rule) {
          h_rule = min(h_rule, rule.h);
        }
      }
      if (h_rule < 1e98) {
        //m_Height[id_node] = h_rule;
      }
    }
  }


  QVector<double> h(m_Grid->GetNumberOfPoints(), 0);
  QVector<double> h_opt(m_Grid->GetNumberOfPoints(), 0);
  int num_layers = 0;
  double err_sq_max = 1e99;
  m_NumLayers = 0;
  do {
    ++num_layers;
    double err_sq = 0;
    double err_max_iter = 0;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_SurfNode[id_node]) {
        double Dh = m_Height[id_node];
        h[id_node] = m_Height[id_node];
        double total_ratio = Dh_max[id_node]/Dh;
        double stretch = pow(total_ratio, 1.0/num_layers);
        for (int i = 1; i < num_layers; ++i) {
          Dh *= stretch;
          h[id_node] += Dh;
        }
        double err;
        if (fabs(stretch) > fabs(m_DesiredStretching)) {
          err = fabs(stretch - m_DesiredStretching)/m_DesiredStretching;
        } else {
          err = fabs(stretch - m_DesiredStretching)/stretch;
        }
        err_max_iter = max(err_max_iter, err);
        err_sq += err*err;
      }
    }
    if (err_sq < err_sq_max) {
      m_NumLayers = num_layers;
      h_opt = h;
      err_sq_max = err_sq;
    }
  } while (num_layers < 200);
  m_Height = h_opt;
}

void GridSmoother::computeHeights()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");

  {
    QString blayer_txt = GuiMainWindow::pointer()->getXmlSection("blayer/global");
    QTextStream s(&blayer_txt);
    if (!s.atEnd()) s >> m_AbsoluteHeight;
    if (!s.atEnd()) s >> m_RelativeHeight;
    if (!s.atEnd()) s >> m_Blending;
    if (!s.atEnd()) s >> m_DesiredStretching;
    if (!s.atEnd()) s >> m_FarRatio;
    if (!s.atEnd()) s >> m_NumLayers;
    if (!s.atEnd()) s >> m_NumHeightRelaxations;
    if (!s.atEnd()) s >> m_NumNormalRelaxations;
  }
  computeDesiredHeights();
  {
    QString blayer_txt = "";
    QTextStream s(&blayer_txt);
    s << m_AbsoluteHeight << " ";
    s << m_RelativeHeight << " ";
    s << m_Blending << " ";
    s << m_DesiredStretching << " ";
    s << m_FarRatio << " ";
    s << m_NumLayers << " ";
    s << m_NumHeightRelaxations << " ";
    s << m_NumNormalRelaxations << " ";
    GuiMainWindow::pointer()->setXmlSection("blayer/global", blayer_txt);
  }

  // compute height limiters
  QVector<double> height_limit(m_Grid->GetNumberOfPoints(), 1e99);

  // gaps
  QList<vtkIdType> search_nodes;
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    bool append_node = m_SurfNode[id_node];
    if (!append_node) {
      for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
        vtkIdType id_cell = m_Part.n2cGG(id_node, i);
        if (isSurface(id_cell, m_Grid)) {
          if (!m_LayerAdjacentBoundaryCodes.contains(bc->GetValue(id_cell))) {
            append_node = true;
            break;
          }
        }
      }
    }
    if (append_node) {
      search_nodes.append(id_node);
    }
  }

  QVector<vec3_t> points(search_nodes.size());
  for (int i = 0; i < search_nodes.size(); ++i) {
    m_Grid->GetPoint(search_nodes[i], points[i].data());
  }
  PointFinder pfind;
  pfind.setPoints(points);

  for (vtkIdType id_node1 = 0; id_node1 < m_Grid->GetNumberOfPoints(); ++id_node1) {
    if (m_SurfNode[id_node1]) {
      vec3_t x1;
      m_Grid->GetPoint(id_node1, x1.data());
      QVector<int> close_points;
      pfind.getClosePoints(x1, close_points, 2*m_Height[id_node1]/m_MaxHeightInGaps);
      foreach (int i, close_points) {
        vtkIdType id_node2 = search_nodes[i];
        if (id_node1 != id_node2) {
          if (m_GeoNormal[id_node1]*m_GeoNormal[id_node2] < 0)  {
            vec3_t x2;
            m_Grid->GetPoint(id_node2, x2.data());
            vec3_t Dx = x2 - x1;
            double a = Dx*m_GeoNormal[id_node1];
            if (a > 0) {
              double c = Dx.abs();
              double b = sqrt(sqr(c) - sqr(a));
              double cos_beta = acos(fabs(m_GeoNormal[id_node1]*m_GeoNormal[id_node2]));
              double e = b/cos_beta;
              double d = sqrt(sqr(e) - sqr(b));
              double alpha = 180.0/M_PI*acos(a/c);
              if (alpha < m_RadarAngle) {
                height_limit[id_node1] = min(height_limit[id_node1], m_MaxHeightInGaps*(a+d));
              }
            }
          }
        }
      }
    }
  }

  // "neighbour" gaps
  for (vtkIdType id_node1 = 0; id_node1 < m_Grid->GetNumberOfPoints(); ++id_node1) {
    if (m_SurfNode[id_node1]) {
      vec3_t x1;
      m_Grid->GetPoint(id_node1, x1.data());
      for (int i = 0; i < m_Part.n2nGSize(id_node1); ++i) {
        vtkIdType id_node2 = m_Part.n2nGG(id_node1, i);
        if (m_SurfNode[id_node2]) {
          vec3_t x2;
          m_Grid->GetPoint(id_node2, x2.data());
          vec3_t v = x2 - x1;
          double L = v.abs();
          v.normalise();
          vec3_t n_plane = v.cross(m_NodeNormal[id_node1]);
          n_plane.normalise(); //??
          vec3_t n1_2d = m_NodeNormal[id_node1] - (m_NodeNormal[id_node1]*n_plane)*n_plane;
          vec3_t n2_2d = m_NodeNormal[id_node2] - (m_NodeNormal[id_node2]*n_plane)*n_plane;
          double s1 = n1_2d*v;
          double s2 = n2_2d*v;
          double m1 = sign1(s1)*sqrt(1 - sqr(s1))/max(1e-3, fabs(s1));
          double m2 = sign1(s2)*sqrt(1 - sqr(s2))/max(1e-3, fabs(s2));
          if (fabs(m1) < 100 && fabs(m2) < 100) {
            double C  = sign1(m2 - m1)*max(1e-3, fabs(m2 - m1));
            double h  = m1*m2/C;
            if (h < 0) h = 1000;
            height_limit[id_node1] = min(height_limit[id_node1], 2*m_MaxHeightInGaps*h*L);
          }
        }
      }
    }
  }

  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    m_Height[id_node] = min(height_limit[id_node], m_Height[id_node]);
  }

  // fifth pass (smoothing)
  QVector<int> num_bcs(m_Grid->GetNumberOfPoints(),0);
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_SurfNode[id_node]) {
      for (int i = 0; i < m_Part.n2bcGSize(id_node); ++i) {
        int bc = m_Part.n2bcG(id_node, i);
        if (m_BoundaryCodes.contains(bc)) {
          ++num_bcs[id_node];
        }
      }
    }
  }
  for (int iter = 0; iter < m_NumHeightRelaxations; ++iter) {
    QVector<double> h_new = m_Height;
    for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
      if (m_SurfNode[id_node]) {
        QList<vtkIdType> snap_points;
        for (int i = 0; i < m_Part.n2nGSize(id_node); ++i) {
          vtkIdType id_neigh = m_Part.n2nGG(id_node,i);
          if (num_bcs[id_node] <= num_bcs[id_neigh]) {
            if (!m_SurfNode[id_neigh]) {
              EG_BUG;
            }
            snap_points.append(id_neigh);
          }
        }
        if (snap_points.size() > 0) {
          h_new[id_node] = 0;
          foreach (vtkIdType id_snap, snap_points) {
            h_new[id_node] += m_Height[id_snap];
          }
          h_new[id_node] /= snap_points.size();
        } else {
          h_new[id_node] = m_Height[id_node];
        }
        //h_new[id_node] = min(h_new[id_node], height_limit[id_node]);
      }
    }
    m_Height = h_new;
  }
}

void GridSmoother::correctNormalVectors()
{
  EG_VTKDCC(vtkIntArray, bc, m_Grid, "cell_code");
  for (vtkIdType id_node = 0; id_node < m_Grid->GetNumberOfPoints(); ++id_node) {
    if (m_NodeNormal[id_node].abs() > 0.1) {
      for (int iter = 0; iter < m_NumBoundaryCorrections; ++iter) {
        for (int i = 0; i < m_Part.n2cGSize(id_node); ++i) {
          vtkIdType id_cell = m_Part.n2cGG(id_node,i);
          if (isSurface(id_cell, m_Grid)) {
            if (!m_BoundaryCodes.contains(bc->GetValue(id_cell))) {
              vec3_t v = m_NodeNormal[id_node];
              vec3_t n = GeometryTools::cellNormal(m_Grid, id_cell);
              n.normalise();
              v -= (n*m_NodeNormal[id_node])*n;
              v.normalise();
              m_NodeNormal[id_node] = v;
            }
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
      m_NodeMarked[pts[0]] = false;
      m_NodeMarked[pts[1]] = false;
      m_NodeMarked[pts[2]] = false;
    }
  }
}

void GridSmoother::simpleNodeMovement(int i_nodes)
{
  vtkIdType id_node = m_Part.globalNode(i_nodes);
  vtkIdType id_foot = m_IdFoot[id_node];
  vec3_t x_new(0,0,0), x_node;
  m_Grid->GetPoint(id_node, x_node.data());

  vec3_t x_surf(0,0,0);
  if (m_SurfNode[id_node]) {
    int N = 0;
    for (int i = 0; i < m_Part.n2nLSize(i_nodes); ++i) {
      vtkIdType id_neigh = m_Part.n2nLG(i_nodes,i);
      if (m_SurfNode[id_neigh]) {
        vec3_t x_neigh;
        m_Grid->GetPoint(id_neigh, x_neigh.data());
        x_surf += x_neigh;
        ++N;
      }
    }
    if (N == 0) {
      EG_BUG;
    } else {
      x_surf *= 1.0/N;
    }
  }

  if (id_foot != -1) {
    vec3_t x_foot, x_node;
    m_Grid->GetPoint(id_foot, x_foot.data());
    //double H = m_Blending*m_AbsoluteHeight + (1.0-m_Blending)*m_RelativeHeight*m_Height[id_foot];
    double H = m_Height[id_foot];
    x_new = x_foot + H*m_NodeNormal[id_foot];
  } else {
    if (m_SurfNode[id_node]) {
      x_new = x_surf;
    } else {
      if (m_Part.n2nLSize(i_nodes) == 0) {
        EG_BUG;
      }
      for (int i = 0; i < m_Part.n2nLSize(i_nodes); ++i) {
        vtkIdType id_neigh = m_Part.n2nLG(i_nodes,i);
        vec3_t x_neigh;
        m_Grid->GetPoint(id_neigh, x_neigh.data());
        x_new += x_neigh;
      }
      x_new *= 1.0/m_Part.n2nLSize(i_nodes);
    }
  }
  vec3_t Dx = x_new - x_node;
  correctDx(i_nodes, Dx);
  moveNode(i_nodes, Dx);
}

void GridSmoother::operate()
{
  if (m_FirstCall) {
    markNodes();
    computeNormals();
    computeFeet();
    computeHeights();
    m_FirstCall = false;
  }
  l2g_t nodes = m_Part.getNodes();
  for (int i_nodes = 0; i_nodes < nodes.size(); ++i_nodes) {
    if (m_NodeMarked[nodes[i_nodes]]) {
      simpleNodeMovement(i_nodes);
    }
  }
}

void GridSmoother::writeDebugFile(QString file_name)
{
  QVector<vtkIdType> bcells;  
  QSet<int> bcs = getAllBoundaryCodes(m_Grid);
  getSurfaceCells(bcs, bcells, m_Grid);
  MeshPartition bpart(m_Grid);
  bpart.setCells(bcells);
  file_name = GuiMainWindow::pointer()->getCwd() + "/" + file_name + ".vtk";
  QFile file(file_name);
  file.open(QIODevice::WriteOnly);
  QTextStream f(&file);
  f << "# vtk DataFile Version 2.0\n";
  f << "m_NodeNormal\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";
  f << "POINTS " << bpart.getNumberOfNodes() << " float\n";
  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vec3_t x;
    vtkIdType id_node = bpart.globalNode(i);
    m_Grid->GetPoint(id_node, x.data());
    f << x[0] << " " << x[1] << " " << x[2] << "\n";
  }
  f << "CELLS " << bpart.getNumberOfCells() << " ";
  int N = 0;
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    N += 1 + N_pts;
  }
  f << N << "\n";
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType N_pts, *pts;
    m_Grid->GetCellPoints(id_cell, N_pts, pts);
    QVector<int> nds(N_pts);
    f << N_pts;
    for (int j = 0; j < N_pts; ++j) {
      f << " " << bpart.localNode(pts[j]);
    }
    f << "\n";
  }

  f << "CELL_TYPES " << bpart.getNumberOfCells() << "\n";
  for (int i = 0; i < bpart.getNumberOfCells(); ++ i) {
    vtkIdType id_cell = bpart.globalCell(i);
    vtkIdType type_cell = m_Grid->GetCellType(id_cell);
    f << type_cell << "\n";
  }
  f << "POINT_DATA " << bpart.getNumberOfNodes() << "\n";
  f << "VECTORS N float\n";

  for (int i = 0; i < bpart.getNumberOfNodes(); ++i) {
    vtkIdType id_node = bpart.globalNode(i);
    f << m_NodeNormal[id_node][0] << " " << m_NodeNormal[id_node][1] << " " << m_NodeNormal[id_node][2] << "\n";
  }
}


